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
#   iteration_results.rds   – list of per-iteration DEA data frames
#   core_markers.csv        – aggregated marker consistency table
#   core_markers.rds        – same data as an R data frame
#   plots/                  – dot plots and heatmaps summarising core markers

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(jsonlite)
  library(ggplot2)
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

set.seed(opt$seed)

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

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

#' Score cells with module scoring and assign labels
assign_labels <- function(sobj, marker_dict, assay) {
  DefaultAssay(sobj) <- assay

  for (ct in names(marker_dict)) {
    genes <- marker_dict[[ct]]
    sobj  <- AddModuleScore(
      sobj,
      features  = list(genes),
      name      = paste0(ct, "_score"),
      assay     = assay,
      seed      = NULL
    )
  }

  score_cols <- make.names(paste0(names(marker_dict), "_score1"))
  score_mat  <- sobj@meta.data[, score_cols, drop = FALSE]
  colnames(score_mat) <- names(marker_dict)

  # Assign each cell to the cell type with the highest module score
  sobj$SCbootstrap_label <- apply(score_mat, 1, function(row) {
    names(marker_dict)[which.max(row)]
  })

  sobj
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

for (i in seq_len(opt$iterations)) {
  message(sprintf("[SCbootstrap] Iteration %d / %d", i, opt$iterations))

  # 5a. Subsample
  sub <- subsample_seurat(seurat_obj, opt$subsample)

  # 5b. Label cells via module scoring
  sub <- assign_labels(sub, marker_dict, opt$assay)

  label_tab <- table(sub$SCbootstrap_label)
  message("  Label distribution: ",
          paste(names(label_tab), label_tab, sep = "=", collapse = "  "))

  # 5c. DEA
  de_result <- run_dea(sub, opt$assay, opt$min_cells, opt$logfc_threshold)

  if (!is.null(de_result)) {
    de_result$iteration <- i
    iteration_results[[i]] <- de_result
  }
}

# Remove empty iterations
iteration_results <- Filter(Negate(is.null), iteration_results)

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
    pct_iter        = round(n() / opt$iterations * 100, 1),
    mean_log2FC     = round(mean(avg_log2FC,  na.rm = TRUE), 4),
    mean_pct1       = round(mean(pct.1,       na.rm = TRUE), 4),
    mean_pct2       = round(mean(pct.2,       na.rm = TRUE), 4),
    mean_p_val_adj  = round(mean(p_val_adj,   na.rm = TRUE), 6),
    .groups = "drop"
  ) %>%
  arrange(cell_type, desc(pct_iter), desc(mean_log2FC))

message("[SCbootstrap] Core markers identified: ", nrow(core_markers), " (gene x cell_type pairs)")

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

# ---------------------------------------------------------------------------
# 8. Plots
# ---------------------------------------------------------------------------

message("[SCbootstrap] Generating summary plots ...")

# 8a. Consistency histogram per cell type
p_hist <- ggplot(core_markers, aes(x = pct_iter, fill = cell_type)) +
  geom_histogram(binwidth = 5, colour = "white", linewidth = 0.2) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(
    title    = "Marker detection consistency across iterations",
    subtitle = sprintf("%d iterations, %.0f%% subsample", opt$iterations,
                       opt$subsample * 100),
    x        = "% iterations detected",
    y        = "Number of genes"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(opt$output_dir, "plots", "consistency_histogram.pdf"),
       p_hist, width = 8, height = 5)

# 8b. Top-N core markers dot plot (top 20 per cell type by consistency)
top_markers <- core_markers %>%
  group_by(cell_type) %>%
  slice_max(order_by = pct_iter, n = 20, with_ties = FALSE) %>%
  ungroup()

p_dot <- ggplot(top_markers,
                aes(x = cell_type, y = reorder(gene, pct_iter),
                    colour = mean_log2FC, size = pct_iter)) +
  geom_point() +
  scale_colour_gradient(low = "grey80", high = "firebrick3",
                        name = "Mean\nlog2FC") +
  scale_size_continuous(range = c(1, 6), name = "% iters\ndetected") +
  labs(
    title = "Top 20 core markers per cell type",
    x     = NULL,
    y     = "Gene"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

n_types <- length(unique(top_markers$cell_type))
n_genes <- length(unique(top_markers$gene))

ggsave(file.path(opt$output_dir, "plots", "top_core_markers_dotplot.pdf"),
       p_dot,
       width  = max(6, n_types * 1.5),
       height = max(5, n_genes * 0.25 + 2))

message("[SCbootstrap] Done.")
message("  iteration_results  -> ", file.path(opt$output_dir, "iteration_results.rds"))
message("  core_markers (CSV) -> ", file.path(opt$output_dir, "core_markers.csv"))
message("  core_markers (RDS) -> ", file.path(opt$output_dir, "core_markers.rds"))
message("  plots/             -> ", file.path(opt$output_dir, "plots"))
