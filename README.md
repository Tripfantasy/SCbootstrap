# SCbootstrap

Iterative subsampling of scRNA-seq data for cell state distillation via bootstrapped
differential expression analysis (DEA).

## Overview

`SCbootstrap.R` takes a Seurat `.RDS` file and a marker dictionary, then:

1. **Subsamples** the Seurat object *n* times (each iteration draws a random fraction of cells).
2. **Labels** cells in each subsample using module scoring against the provided marker dictionary.
3. **Runs lightweight DEA** (Wilcoxon rank-sum test via `FindMarkers`) per cell label per iteration.
4. **Aggregates** all DEA results across iterations to identify **core markers** – genes
   that are consistently differentially expressed for each cell type.

## Requirements

```r
install.packages(c("optparse", "dplyr", "jsonlite", "ggplot2"))
# Seurat v4 or v5
install.packages("Seurat")
```

## Usage

```bash
Rscript SCbootstrap.R \
  --input        path/to/input.rds   \
  --markers      path/to/markers.json \
  --iterations   50                  \
  --subsample    0.8                 \
  --output_dir   results/
```

### All options

| Flag | Default | Description |
|---|---|---|
| `--input` | *(required)* | Path to input Seurat `.RDS` file |
| `--markers` | *(required)* | Marker dictionary: `.json` or `.csv` (see below) |
| `--iterations` | `50` | Number of bootstrap iterations |
| `--subsample` | `0.8` | Fraction of cells sampled per iteration |
| `--output_dir` | `SCbootstrap_results` | Output directory (created if absent) |
| `--assay` | `RNA` | Seurat assay to use |
| `--min_cells` | `10` | Min cells per label to run DEA |
| `--logfc_threshold` | `0.25` | Log-fold-change threshold for `FindMarkers` |
| `--seed` | `42` | Random seed |

## Examples

### Minimal run

The only required arguments are `--input` and `--markers`. All other flags use
their defaults (50 iterations, 80 % subsample, output to `SCbootstrap_results/`).

```bash
Rscript SCbootstrap.R \
  --input   data/pbmc3k.rds \
  --markers markers/pbmc_markers.json
```

---

### Quick sanity check (5 iterations, small subsample)

Run a fast sanity check before committing to a full analysis. Using a small
subsample (`0.3`) and few iterations (`5`) completes in seconds on a laptop.

```bash
Rscript SCbootstrap.R \
  --input        data/pbmc3k.rds \
  --markers      markers/pbmc_markers.json \
  --iterations   5 \
  --subsample    0.3 \
  --output_dir   sanity_check/
```

---

### Full run with custom parameters

Override every default to increase statistical power and tighten thresholds:

```bash
Rscript SCbootstrap.R \
  --input            data/pbmc3k.rds \
  --markers          markers/pbmc_markers.json \
  --iterations       100 \
  --subsample        0.7 \
  --assay            RNA \
  --min_cells        20 \
  --logfc_threshold  0.5 \
  --output_dir       results/pbmc_100iter/
```

---

### Reproducible run (fixed seed)

Set `--seed` to an explicit value so the exact same cells are drawn in every
iteration, making results fully reproducible:

```bash
Rscript SCbootstrap.R \
  --input      data/pbmc3k.rds \
  --markers    markers/pbmc_markers.json \
  --seed       123 \
  --output_dir results/pbmc_seed123/
```

---

### Using a CSV marker file

If you prefer a CSV over JSON, provide the path to the `.csv` file.
The file must have two columns (no header): `CellType`, `Gene`.

```bash
Rscript SCbootstrap.R \
  --input      data/pbmc3k.rds \
  --markers    markers/pbmc_markers.csv \
  --output_dir results/pbmc_csv/
```

---

### Reading results in R

After the script finishes, load the outputs directly in R for further analysis:

```r
# Core markers table (aggregated across all iterations)
core <- read.csv("results/pbmc_100iter/core_markers.csv")

# View top markers per cell type (>= 80 % iteration detection)
subset(core, pct_iter >= 80)

# Per-iteration DEA data frames (list, one element per iteration)
iter_results <- readRDS("results/pbmc_100iter/iteration_results.rds")
head(iter_results[[1]])   # inspect the first iteration
```

---

## Marker dictionary formats

**JSON** (recommended)
```json
{
  "T_cell":    ["CD3D", "CD3E", "CD8A"],
  "B_cell":    ["CD19", "MS4A1", "CD79A"],
  "Monocyte":  ["CD14", "LYZ", "CST3"]
}
```

**CSV** (no header, two columns: CellType, Gene)
```
T_cell,CD3D
T_cell,CD3E
B_cell,CD19
```

## Outputs

All files are written to `--output_dir`:

| File | Description |
|---|---|
| `iteration_results.rds` | List of per-iteration DEA data frames |
| `core_markers.csv` | Aggregated marker consistency table (one row per gene × cell type) |
| `core_markers.rds` | Same table as an R data frame |
| `plots/consistency_histogram.pdf` | Detection frequency histogram per cell type |
| `plots/top_core_markers_dotplot.pdf` | Dot plot of top-20 core markers per cell type |

### `core_markers.csv` columns

| Column | Description |
|---|---|
| `cell_type` | Cell type label |
| `gene` | Gene name |
| `n_iter_detected` | Number of iterations in which this gene was a significant DE marker |
| `pct_iter` | Percentage of iterations detected |
| `mean_log2FC` | Mean log₂ fold-change across detected iterations |
| `mean_pct1` | Mean fraction expressing in the labelled cell type |
| `mean_pct2` | Mean fraction expressing in background cells |
| `mean_p_val_adj` | Mean adjusted p-value across iterations |
