#!/usr/bin/env Rscript
# SCbootstrap.R
# Iterative subsampling of scRNA-seq data for cell state distillation.
#
# Usage:
#   Rscript SCbootstrap.R \
#     --input        path/to/input.rds        \
#     --markers      path/to/markers.json      \
#     --iterations   50                        \
#     --subsample    0.8                        \
#     --output_dir   results/                  \
#     [--assay       RNA]                      \
#     [--min_cells   10]                       \
#     [--logfc_threshold 0.25]                 \
#     [--seed        42]
#
# Input:
#   --input        : Path to a Seurat object saved as .RDS
#   --markers      : Path to a JSON or CSV file containing the marker dictionary.
#                    JSON format: { "CellType1": ["GeneA","GeneB"], "CellType2": [...] }
#                    CSV format : two columns, no header required – CellType, Gene
#   --iterations   : Number of bootstrap iterations (default: 50)
#   --subsample    : Fraction of cells to sample per iteration (default: 0.8)
#   --output_dir   : Directory to write results (created if absent)
#   --assay        : Seurat assay to use (default: "RNA")
#   --min_cells    : Minimum number of cells required per label to run DEA (default: 10)
#   --logfc_threshold : Log-fold-change threshold for FindMarkers (default: 0.25)
#   --seed         : Random seed for reproducibility (default: 42)
#
# Outputs (written to output_dir):
#   iteration_results.rds          – list of per-iteration DEA data frames
#   core_markers.csv               – aggregated marker consistency table
#   core_markers.rds               – same data as an R data frame
#   core_gene_metrics.csv          – gene-level stability/contribution metrics
#   core_gene_sets.rds             – core genes per cell type
#   iteration_celltype_metrics.csv – per-iteration label/distribution diagnostics
#   interactive_report.html        – interactive marker contribution report

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(jsonlite)
  library(plotly)
  library(DT)
  library(htmltools)
})

# ---------------------------------------------------------------------------
# 1. Argument parsing
# ---------------------------------------------------------------------------

option_list <- list(
  make_option("--input",           type = "character", default = NULL,
              help = "Path to input Seurat .RDS file [required]"),
  make_option("--markers",         type = "character", default = NULL,
              help = "Path to marker dictionary (JSON or CSV) [required]"),
  make_option("--iterations",      type = "integer",   default = 50L,
              help = "Number of bootstrap iterations [default: %default]"),
  make_option("--subsample",       type = "double",    default = 0.8,
              help = "Fraction of cells to subsample per iteration [default: %default]"),
  make_option("--output_dir",      type = "character", default = "SCbootstrap_results",
              help = "Output directory [default: %default]"),
  make_option("--assay",           type = "character", default = "RNA",
              help = "Seurat assay to use [default: %default]"),
  make_option("--min_cells",       type = "integer",   default = 10L,
              help = "Minimum cells per label for DEA [default: %default]"),
  make_option("--logfc_threshold", type = "double",    default = 0.25,
              help = "Log-fold-change threshold for FindMarkers [default: %default]"),
  make_option("--core_pct_threshold", type = "double", default = 70,
              help = "Min %% detection for core gene calling [default: %default]"),
  make_option("--core_cv_threshold", type = "double", default = 1.0,
              help = "Max coefficient of variation for core gene calling [default: %default]"),
  make_option("--seed",            type = "integer",   default = 42L,
              help = "Random seed [default: %default]")
)

parse_args_safe <- function(option_list) {
  parser <- OptionParser(option_list = option_list,
                         description = paste(
                           "SCbootstrap: iterative subsampling of scRNA-seq data",
                           "for cell state distillation."
                         ))
  tryCatch(
    parse_args(parser),
    error = function(e) {
      message("Argument error: ", conditionMessage(e))
      print_help(parser)
      quit(status = 1)
    }
  )
}

opt <- parse_args_safe(option_list)

if (is.null(opt$input))   stop("--input is required.")
if (is.null(opt$markers)) stop("--markers is required.")
if (!file.exists(opt$input))   stop("Input file not found: ", opt$input)
if (!file.exists(opt$markers)) stop("Markers file not found: ", opt$markers)
if (opt$subsample <= 0 || opt$subsample > 1)
  stop("--subsample must be in (0, 1].")
if (opt$iterations < 1)
  stop("--iterations must be >= 1.")
if (opt$core_pct_threshold < 0 || opt$core_pct_threshold > 100)
  stop("--core_pct_threshold must be in [0, 100].")
if (opt$core_cv_threshold < 0)
  stop("--core_cv_threshold must be >= 0.")

set.seed(opt$seed)

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------------------------

message("[SCbootstrap] Loading Seurat object from: ", opt$input)
seurat_obj <- readRDS(opt$input)

if (!inherits(seurat_obj, "Seurat")) {
  stop("The .RDS file does not contain a Seurat object.")
}

# Ensure the requested assay is present
if (!(opt$assay %in% names(seurat_obj@assays))) {
  stop("Assay '", opt$assay, "' not found in the Seurat object. ",
       "Available assays: ", paste(names(seurat_obj@assays), collapse = ", "))
}

DefaultAssay(seurat_obj) <- opt$assay
message("[SCbootstrap] Cells: ", ncol(seurat_obj),
        "  |  Genes: ", nrow(seurat_obj))

# ---------------------------------------------------------------------------
# 3. Load marker dictionary
# ---------------------------------------------------------------------------

load_markers <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "json") {
    raw <- jsonlite::read_json(path, simplifyVector = TRUE)
    # Accepts list of character vectors
    lapply(raw, as.character)
  } else if (ext %in% c("csv", "tsv")) {
    sep <- if (ext == "tsv") "\t" else ","
    df  <- read.csv(path, header = FALSE, sep = sep, stringsAsFactors = FALSE,
                    col.names = c("CellType", "Gene"))
    split(df$Gene, df$CellType)
  } else {
    stop("Unsupported marker file format: '", ext,
         "'. Please provide a .json, .csv, or .tsv file.")
  }
}

message("[SCbootstrap] Loading marker dictionary from: ", opt$markers)
marker_dict <- load_markers(opt$markers)

# Filter genes in marker dict to those present in the object
all_genes <- rownames(seurat_obj)
marker_dict <- lapply(marker_dict, function(genes) {
  present <- genes[genes %in% all_genes]
  if (length(present) == 0) NULL else present
})
marker_dict <- Filter(Negate(is.null), marker_dict)

if (length(marker_dict) == 0) {
  stop("No marker genes from the dictionary were found in the Seurat object.")
}

cell_types <- names(marker_dict)
message("[SCbootstrap] Cell types in dictionary: ", paste(cell_types, collapse = ", "))

# ---------------------------------------------------------------------------
# 4. Helper functions
# ---------------------------------------------------------------------------

#' Subsample a Seurat object
subsample_seurat <- function(sobj, frac) {
  n_cells  <- ncol(sobj)
  n_sample <- max(1L, floor(n_cells * frac))
  sampled  <- sample(colnames(sobj), n_sample, replace = FALSE)
  sobj[, sampled]
}

#' Get SCbootstrap module score columns
get_score_columns <- function(marker_dict) {
  setNames(
   paste0(make.names(names(marker_dict)), "_score1"),
   names(marker_dict)
  )
}

#' Score cells with module scoring and assign labels
assign_labels <- function(sobj, marker_dict, assay) {
  DefaultAssay(sobj) <- assay

  for (ct in names(marker_dict)) {
   genes <- marker_dict[[ct]]
   sobj  <- AddModuleScore(
     sobj,
     features  = list(genes),
     name      = paste0(make.names(ct), "_score"),
     assay     = assay,
     seed      = NULL
   )
  }

  score_cols <- get_score_columns(marker_dict)
  score_mat  <- sobj@meta.data[, unname(score_cols), drop = FALSE]
  colnames(score_mat) <- names(marker_dict)

  # Assign each cell to the cell type with the highest module score
  sobj$SCbootstrap_label <- apply(score_mat, 1, function(row) {
   names(marker_dict)[which.max(row)]
  })

  sobj
}

#' Summarise marker contribution and iteration distributions
build_iteration_diagnostics <- function(sobj, marker_dict, iteration, min_cells) {
  score_cols <- get_score_columns(marker_dict)
  score_mat  <- sobj@meta.data[, unname(score_cols), drop = FALSE]
  colnames(score_mat) <- names(marker_dict)
  score_df <- cbind(
    data.frame(cell_id = rownames(score_mat),
               assigned_label = sobj$SCbootstrap_label,
               stringsAsFactors = FALSE),
    as.data.frame(score_mat, stringsAsFactors = FALSE)
  )

  contrib_rows <- lapply(names(marker_dict), function(ct) {
    vals <- score_df[[ct]]
    in_ct <- score_df$assigned_label == ct
    data.frame(
      iteration = iteration,
      cell_type = ct,
      marker_set = ct,
      mean_score_all = mean(vals, na.rm = TRUE),
      mean_score_assigned = if (sum(in_ct) > 0) mean(vals[in_ct], na.rm = TRUE) else NA_real_,
      median_score_assigned = if (sum(in_ct) > 0) median(vals[in_ct], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  contribution_metrics <- do.call(rbind, contrib_rows)

  labels <- names(marker_dict)
  label_counts <- vapply(labels, function(lbl) sum(sobj$SCbootstrap_label == lbl), integer(1))
  iter_metrics <- data.frame(
    iteration = iteration,
    cell_type = labels,
    n_cells = as.integer(label_counts),
    stringsAsFactors = FALSE
  )
  iter_metrics$de_tested <- vapply(seq_len(nrow(iter_metrics)), function(i) {
    n_fg <- iter_metrics$n_cells[i]
    n_bg <- sum(iter_metrics$n_cells) - n_fg
    n_fg >= min_cells && n_bg >= min_cells
  }, logical(1))

  list(
    contribution_metrics = contribution_metrics,
    iteration_metrics = iter_metrics
  )
}

#' Run lightweight DEA (Wilcoxon) for one labelled Seurat object
run_dea <- function(sobj, assay, min_cells, logfc_threshold) {
  Idents(sobj) <- "SCbootstrap_label"
  labels       <- levels(factor(sobj$SCbootstrap_label))

  results <- list()
  for (lbl in labels) {
    n_fg <- sum(sobj$SCbootstrap_label == lbl)
    n_bg <- sum(sobj$SCbootstrap_label != lbl)

    if (n_fg < min_cells || n_bg < min_cells) {
      message("  Skipping label '", lbl, "' (fg=", n_fg, ", bg=", n_bg, ")")
      next
    }

    #' Compute robust core-gene metrics from aggregate DEA output
    identify_core_genes <- function(core_markers, all_de, core_pct_threshold, core_cv_threshold, min_log2fc) {
      stability <- all_de %>%
        group_by(cell_type, gene) %>%
        summarise(
          sd_log2FC = sd(avg_log2FC, na.rm = TRUE),
          .groups = "drop"
        )

      gene_metrics <- core_markers %>%
        left_join(stability, by = c("cell_type", "gene")) %>%
        mutate(
          cv_log2FC = ifelse(abs(mean_log2FC) > 1e-12, sd_log2FC / abs(mean_log2FC), NA_real_),
          cv_log2FC = ifelse(is.nan(cv_log2FC) | is.infinite(cv_log2FC), NA_real_, cv_log2FC),
          contribution_score = round((pct_iter * pmax(mean_log2FC, 0)) / (1 + coalesce(cv_log2FC, 1)), 4),
          is_core_gene = pct_iter >= core_pct_threshold &
            mean_log2FC >= min_log2fc &
            coalesce(cv_log2FC, 0) <= core_cv_threshold
        ) %>%
        arrange(cell_type, desc(is_core_gene), desc(contribution_score), desc(mean_log2FC))

      cell_types <- unique(gene_metrics$cell_type)
      core_gene_sets <- setNames(lapply(cell_types, function(ct) {
        gene_metrics$gene[gene_metrics$cell_type == ct & gene_metrics$is_core_gene]
      }), cell_types)

      celltype_summary <- gene_metrics %>%
        group_by(cell_type) %>%
        summarise(
          n_genes = n(),
          n_core_genes = sum(is_core_gene),
          median_pct_iter = round(median(pct_iter, na.rm = TRUE), 2),
          median_contribution_score = round(median(contribution_score, na.rm = TRUE), 4),
          .groups = "drop"
        )

      list(
        gene_metrics = gene_metrics,
        core_gene_sets = core_gene_sets,
        celltype_summary = celltype_summary
      )
    }

    #' Generate interactive HTML report
    generate_interactive_report <- function(core_markers, core_gene_metrics, celltype_summary,
                                            iteration_metrics, contribution_metrics, output_file) {
      top_markers <- core_gene_metrics %>%
        group_by(cell_type) %>%
        slice_max(order_by = contribution_score, n = 25, with_ties = FALSE) %>%
        ungroup()

      p_contrib <- plot_ly(
        data = top_markers,
        x = ~cell_type,
        y = ~gene,
        type = "scatter",
        mode = "markers",
        color = ~mean_log2FC,
        colors = "Reds",
        marker = list(sizemode = "diameter"),
        size = ~pmax(pct_iter, 1),
        text = ~paste0(
          "Cell type: ", cell_type,
          "<br>Gene: ", gene,
          "<br>% detected: ", pct_iter,
          "<br>Mean log2FC: ", mean_log2FC,
          "<br>Contribution score: ", contribution_score,
          "<br>Core gene: ", is_core_gene
        ),
        hoverinfo = "text"
      ) %>%
        layout(
          title = "Marker contribution map (top genes by contribution score)",
          xaxis = list(title = "Cell type"),
          yaxis = list(title = "Gene")
        )

      p_distribution <- plot_ly(
        data = core_markers,
        x = ~cell_type,
        y = ~pct_iter,
        type = "box",
        color = ~cell_type,
        boxpoints = "outliers",
        showlegend = FALSE
      ) %>%
        layout(
          title = "Distribution of marker detection consistency by cell type",
          xaxis = list(title = "Cell type"),
          yaxis = list(title = "% iterations detected")
        )

      contribution_table <- contribution_metrics %>%
        mutate(
          mean_score_all = round(mean_score_all, 4),
          mean_score_assigned = round(mean_score_assigned, 4),
          median_score_assigned = round(median_score_assigned, 4)
        )

      iter_table <- iteration_metrics %>%
        arrange(cell_type, iteration)

      report_tag <- tagList(
        tags$h1("SCbootstrap interactive report"),
        tags$p("Inspect marker contributions, cell-type distribution metrics, and core-gene consistency."),
        tags$h2("Marker contribution map"),
        as.tags(p_contrib),
        tags$h2("Detection consistency distribution"),
        as.tags(p_distribution),
        tags$h2("Cell-type summary"),
        DT::datatable(celltype_summary, options = list(pageLength = 10)),
        tags$h2("Marker contribution metrics by iteration"),
        DT::datatable(contribution_table, options = list(pageLength = 10)),
        tags$h2("Per-iteration cell-type distribution metrics"),
        DT::datatable(iter_table, options = list(pageLength = 10)),
        tags$h2("Gene-level metrics"),
        DT::datatable(
          core_gene_metrics %>%
            select(cell_type, gene, pct_iter, mean_log2FC, cv_log2FC, contribution_score, is_core_gene),
          options = list(pageLength = 15)
        )
      )

      htmltools::save_html(
        html = browsable(report_tag),
        file = output_file
      )
    }

    de <- tryCatch(
      FindMarkers(
        sobj,
        ident.1          = lbl,
        assay            = assay,
        test.use         = "wilcox",
        logfc.threshold  = logfc_threshold,
        min.pct          = 0.1,
        only.pos         = TRUE,
        verbose          = FALSE
      ),
      error = function(e) {
        message("  DEA failed for '", lbl, "': ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(de) && nrow(de) > 0) {
      de$gene      <- rownames(de)
      de$cell_type <- lbl
      results[[lbl]] <- de
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# 5. Main bootstrap loop
# ---------------------------------------------------------------------------

message("[SCbootstrap] Starting ", opt$iterations, " iterations ...")

iteration_results <- vector("list", opt$iterations)
iteration_diagnostics <- vector("list", opt$iterations)
iteration_contributions <- vector("list", opt$iterations)

for (i in seq_len(opt$iterations)) {
  message(sprintf("[SCbootstrap] Iteration %d / %d", i, opt$iterations))

  # 5a. Subsample
  sub <- subsample_seurat(seurat_obj, opt$subsample)

  # 5b. Label cells via module scoring
  sub <- assign_labels(sub, marker_dict, opt$assay)

  label_tab <- table(sub$SCbootstrap_label)
  message("  Label distribution: ",
          paste(names(label_tab), label_tab, sep = "=", collapse = "  "))

  diagnostics <- build_iteration_diagnostics(
    sobj = sub,
    marker_dict = marker_dict,
    iteration = i,
    min_cells = opt$min_cells
  )
  iteration_diagnostics[[i]] <- diagnostics$iteration_metrics
  iteration_contributions[[i]] <- diagnostics$contribution_metrics

  # 5c. DEA
  de_result <- run_dea(sub, opt$assay, opt$min_cells, opt$logfc_threshold)

  if (!is.null(de_result)) {
    de_result$iteration <- i
    iteration_results[[i]] <- de_result
  }
}

# Remove empty iterations
iteration_results <- Filter(Negate(is.null), iteration_results)
iteration_diagnostics <- Filter(Negate(is.null), iteration_diagnostics)
iteration_contributions <- Filter(Negate(is.null), iteration_contributions)

if (length(iteration_results) == 0) {
  stop("No DEA results were produced. Consider lowering --min_cells or --logfc_threshold.")
}

# ---------------------------------------------------------------------------
# 6. Aggregate results – identify core markers
# ---------------------------------------------------------------------------

message("[SCbootstrap] Aggregating results ...")

all_de <- do.call(rbind, iteration_results)
rownames(all_de) <- NULL

# For each (cell_type, gene), compute:
#   - n_iter_detected : number of iterations in which the gene was detected
#   - mean_log2FC     : mean log2 fold-change across those iterations
#   - mean_pct1       : mean fraction of cells expressing in cell type
#   - mean_pct2       : mean fraction of cells expressing in background
#   - mean_p_val_adj  : mean adjusted p-value

core_markers <- all_de %>%
  group_by(cell_type, gene) %>%
  summarise(
    n_iter_detected = n(),
    mean_log2FC     = round(mean(avg_log2FC,  na.rm = TRUE), 4),
    mean_pct1       = round(mean(pct.1,       na.rm = TRUE), 4),
    mean_pct2       = round(mean(pct.2,       na.rm = TRUE), 4),
    mean_p_val_adj  = round(mean(p_val_adj,   na.rm = TRUE), 6),
    .groups = "drop"
  )

celltype_iter_metrics <- do.call(rbind, iteration_diagnostics)
celltype_denominators <- celltype_iter_metrics %>%
  group_by(cell_type) %>%
  summarise(
    n_iter_label_present = sum(n_cells > 0),
    n_iter_de_tested = sum(de_tested),
    mean_cells_per_iteration = round(mean(n_cells), 2),
    sd_cells_per_iteration = round(sd(n_cells), 2),
    .groups = "drop"
  )

core_markers <- core_markers %>%
  left_join(celltype_denominators, by = "cell_type") %>%
  mutate(
    pct_iter = round(ifelse(n_iter_de_tested > 0, n_iter_detected / n_iter_de_tested * 100, 0), 1)
  ) %>%
  arrange(cell_type, desc(pct_iter), desc(mean_log2FC))

message("[SCbootstrap] Core markers identified: ", nrow(core_markers), " (gene x cell_type pairs)")

core_gene_out <- identify_core_genes(
  core_markers = core_markers,
  all_de = all_de,
  core_pct_threshold = opt$core_pct_threshold,
  core_cv_threshold = opt$core_cv_threshold,
  min_log2fc = opt$logfc_threshold
)

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------

message("[SCbootstrap] Writing outputs to: ", opt$output_dir)

saveRDS(iteration_results,
        file.path(opt$output_dir, "iteration_results.rds"))

saveRDS(core_markers,
        file.path(opt$output_dir, "core_markers.rds"))

write.csv(core_markers,
          file.path(opt$output_dir, "core_markers.csv"),
          row.names = FALSE)

write.csv(core_gene_out$gene_metrics,
          file.path(opt$output_dir, "core_gene_metrics.csv"),
          row.names = FALSE)

write.csv(celltype_iter_metrics,
          file.path(opt$output_dir, "iteration_celltype_metrics.csv"),
          row.names = FALSE)

saveRDS(core_gene_out$core_gene_sets,
        file.path(opt$output_dir, "core_gene_sets.rds"))

# ---------------------------------------------------------------------------
# 8. Interactive report
# ---------------------------------------------------------------------------

message("[SCbootstrap] Generating interactive report ...")
contribution_metrics <- do.call(rbind, iteration_contributions)
generate_interactive_report(
  core_markers = core_markers,
  core_gene_metrics = core_gene_out$gene_metrics,
  celltype_summary = core_gene_out$celltype_summary,
  iteration_metrics = celltype_iter_metrics,
  contribution_metrics = contribution_metrics,
  output_file = file.path(opt$output_dir, "interactive_report.html")
)

message("[SCbootstrap] Done.")
message("  iteration_results  -> ", file.path(opt$output_dir, "iteration_results.rds"))
message("  core_markers (CSV) -> ", file.path(opt$output_dir, "core_markers.csv"))
message("  core_markers (RDS) -> ", file.path(opt$output_dir, "core_markers.rds"))
message("  core_gene_metrics  -> ", file.path(opt$output_dir, "core_gene_metrics.csv"))
message("  core_gene_sets     -> ", file.path(opt$output_dir, "core_gene_sets.rds"))
message("  iter_metrics       -> ", file.path(opt$output_dir, "iteration_celltype_metrics.csv"))
message("  interactive report -> ", file.path(opt$output_dir, "interactive_report.html"))
