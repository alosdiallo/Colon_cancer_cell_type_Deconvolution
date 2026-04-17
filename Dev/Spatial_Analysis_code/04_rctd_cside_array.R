#!/usr/bin/env Rscript
# ============================================================
# 04_rctd_cside_array.R
#
# Processes ONE Visium sample with RCTD deconvolution + C-SIDE
# cell-type-specific differential expression.
# Called by submit_rctd_cside_array.sh
#
# Reads from prepped_samples/ (same RDS files used by CellChat)
# to avoid Python/h5ad dependency on compute nodes.
#
# Usage: Rscript 04_rctd_cside_array.R <TASK_ID>
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 04_rctd_cside_array.R <TASK_ID>")
task_id <- as.integer(args[1])

# ============================================================
# HPC LIBRARIES
# ============================================================
.libPaths(c(
  "/dartfs-hpc/rc/home/k/f006fpk/R/x86_64-redhat-linux-gnu-library/4.5",
  "/dartfs-hpc/admin/opt/el8/Rlibs/4.5",
  "/dartfs-hpc/admin/opt/el8/R/4.5.1/lib/R/library"
))

suppressPackageStartupMessages({
  library(spacexr)
  library(Seurat)
  library(Matrix)
})

# ============================================================
# PARALLEL SETTINGS
# The atlas reference is ~14GB. Each RCTD parallel worker loads
# a copy, so max_cores is capped at 4 to stay within 120G.
# (main process ~20GB + 4 workers × ~14GB ≈ 76GB)
# ============================================================
ncores <- 4
cat("Running RCTD with max_cores =", ncores, "\n")

# ============================================================
# PATHS
# ============================================================
base_dir   <- "/dartfs/rc/nosnapshots/V/VaickusL-nb/EDIT_Students/users/alos/colon_work/adata"
prep_dir   <- file.path(base_dir, "prepped_samples")
atlas_path <- file.path(base_dir, "atlas_data.rds")
output_dir <- file.path(base_dir, "rctd_cside_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# FIND SAMPLE FOR THIS TASK
# ============================================================
prep_files <- sort(list.files(prep_dir, pattern = "^prep_.*\\.rds$", full.names = FALSE))
sample_ids <- gsub("^prep_|\\.rds$", "", prep_files)

cat("Found", length(sample_ids), "prepped samples\n")
cat("Sample IDs:", paste(sample_ids, collapse = ", "), "\n\n")

if (task_id >= length(sample_ids)) {
  cat("Task ID", task_id, "exceeds number of samples (", length(sample_ids), "). Exiting.\n")
  quit(status = 0)
}

sid <- sample_ids[task_id + 1]   # +1: R is 1-indexed, SLURM array is 0-indexed
cat("Task ID:", task_id, "| Sample:", sid, "\n")

# Skip if already completed successfully
out_rctd  <- file.path(output_dir, paste0("rctd_", sid, ".rds"))
out_cside <- file.path(output_dir, paste0("cside_", sid, ".rds"))

# Check if C-SIDE result is real or an error placeholder
cside_done <- FALSE
if (file.exists(out_cside)) {
  cside_check <- readRDS(out_cside)
  if (!is.list(cside_check) || !"status" %in% names(cside_check)) {
    cside_done <- TRUE   # Real C-SIDE result, not an error object
  } else {
    cat("Previous C-SIDE was:", cside_check$status, "— will re-run.\n")
    file.remove(out_cside)
  }
  rm(cside_check)
}

if (file.exists(out_rctd) && cside_done) {
  cat("Sample", sid, "already complete (both RCTD and C-SIDE). Skipping.\n")
  quit(status = 0)
}

# ============================================================
# STEP 1: BUILD RCTD REFERENCE FROM ATLAS
# ============================================================
cat("\n=== Step 1: Building Reference ===\n")

if (!file.exists(atlas_path)) stop("Cannot find atlas_data.rds at: ", atlas_path)

atlas <- readRDS(atlas_path)

# Extract raw counts (RCTD requires raw counts, not normalized)
counts_ref <- tryCatch(
  GetAssayData(atlas, layer = "counts"),
  error = function(e) atlas@assays$RNA@counts
)
if (ncol(counts_ref) == 0 || nrow(counts_ref) == 0) {
  counts_ref <- atlas@assays$RNA@counts
}
cat("  Atlas dimensions:", nrow(counts_ref), "genes x", ncol(counts_ref), "cells\n")

# Cell type annotations
cell_types <- as.factor(atlas@meta.data$cl295v11SubFull)
names(cell_types) <- colnames(atlas)
cat("  Cell types:", nlevels(cell_types), "unique types\n")

# Sanitize cell type names: replace / with _ (spacexr chokes on /)
cell_types <- as.factor(gsub("/", "_", as.character(cell_types)))
names(cell_types) <- colnames(counts_ref)

nUMI_ref <- colSums(counts_ref)
reference <- Reference(counts_ref, cell_types, nUMI_ref)
cat("  Reference created successfully\n")

rm(atlas); gc()

# ============================================================
# STEP 2: LOAD VISIUM SAMPLE FROM PREPPED RDS
# ============================================================
cat("\n=== Step 2: Loading prepped sample", sid, "===\n")

prep <- readRDS(file.path(prep_dir, paste0("prep_", sid, ".rds")))

# Extract raw counts from the Seurat object
counts_vis <- tryCatch(
  GetAssayData(prep$seu, layer = "counts", assay = "RNA"),
  error = function(e) GetAssayData(prep$seu, slot = "counts", assay = "RNA")
)
if (!inherits(counts_vis, "dgCMatrix")) {
  counts_vis <- as(counts_vis, "dgCMatrix")
}
cat("  Visium dimensions:", nrow(counts_vis), "genes x", ncol(counts_vis), "spots\n")

# Spatial coordinates
coords <- as.data.frame(prep$coords)
if (ncol(coords) >= 2) {
  colnames(coords)[1:2] <- c("x", "y")
}
# Ensure rownames match count matrix columns
rownames(coords) <- colnames(counts_vis)
cat("  Coordinates loaded:", nrow(coords), "spots\n")

nUMI_vis <- colSums(counts_vis)
spatial_rna <- SpatialRNA(coords, counts_vis, nUMI_vis)
cat("  SpatialRNA created\n")

# ============================================================
# STEP 3: RUN RCTD
# ============================================================
cat("\n=== Step 3: Running RCTD ===\n")

if (file.exists(out_rctd)) {
  cat("  RCTD output already exists. Loading from disk.\n")
  rctd <- readRDS(out_rctd)
} else {
  rctd <- create.RCTD(spatial_rna, reference, max_cores = 1)
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  saveRDS(rctd, out_rctd)
  cat("  Saved:", out_rctd, "\n")
}

cat("  RCTD complete\n")

# Normalize weights to proportions
weights <- rctd@results$weights
weights_norm <- sweep(weights, 1, rowSums(weights), "/")
cat("  Weight matrix:", nrow(weights), "spots x", ncol(weights), "cell types\n")
cat("  Mean proportions (top 10):\n")
print(head(sort(colMeans(weights_norm), decreasing = TRUE), 10))

# ============================================================
# STEP 4: RUN C-SIDE
# ============================================================
cat("\n=== Step 4: Running C-SIDE ===\n")

# Extract histology covariate from prepped metadata
histology <- prep$meta$histology

if (is.null(histology)) {
  cat("  WARNING: No histology column found in sample", sid, ". Skipping C-SIDE.\n")
  saveRDS(list(sample_id = sid, status = "skipped_no_histology"), out_cside)
  quit(status = 0)
}

# Binary covariate: cancer tissue (CT) = 1, non-cancer = 0
region_binary <- as.numeric(histology == "cancer_cancer")
names(region_binary) <- rownames(prep$meta)
n_ct    <- sum(region_binary == 1)
n_nonct <- sum(region_binary == 0)
balance_ratio <- min(n_ct, n_nonct) / max(n_ct, n_nonct)
cat("  Region covariate: CT =", n_ct, "spots, non-CT =", n_nonct, "spots\n")
cat("  Balance ratio:", round(balance_ratio, 3), "\n")

# Check for degenerate case: need spots in both groups
if (n_ct < 10 || n_nonct < 10) {
  cat("  WARNING: Insufficient spots in one region group. Skipping C-SIDE.\n")
  saveRDS(list(sample_id = sid, status = "skipped_insufficient_regions",
               n_ct = n_ct, n_nonct = n_nonct, balance_ratio = balance_ratio),
          out_cside)
  quit(status = 0)
}

# Skip severely imbalanced samples — they will still be included in
# run.CSIDE.replicates for the population-level mets comparison
if (balance_ratio < 0.1) {
  cat("  WARNING: Severe imbalance (ratio =", round(balance_ratio, 3),
      "). Skipping single-sample C-SIDE.\n")
  cat("  This sample will still be included in run.CSIDE.replicates.\n")
  saveRDS(list(sample_id = sid, status = "skipped_imbalanced",
               n_ct = n_ct, n_nonct = n_nonct, balance_ratio = balance_ratio),
          out_cside)
  quit(status = 0)
}

# Let run.CSIDE.single determine eligible cell types internally
# via aggregate_cell_types. Our manual pre-filter was inconsistent
# with the internal filter, causing errors.
cat("  Cell types in RCTD:", ncol(weights_norm), "\n")
cat("  Spots above 5% per type (top 15):\n")
spots_above <- sort(colSums(weights_norm > 0.05), decreasing = TRUE)
print(head(spots_above, 15))

# Align explanatory variable to RCTD barcodes
barcodes_use <- colnames(rctd@spatialRNA@counts)
ev <- region_binary[barcodes_use]

tryCatch({
  rctd_cside <- run.CSIDE.single(
    myRCTD                = rctd,
    explanatory.variable  = ev,
    cell_type_threshold   = 25,
    doublet_mode          = FALSE,   # MUST be FALSE for full-mode RCTD
    fdr                   = 0.05,
    weight_threshold      = 0.75
  )

  cat("  C-SIDE complete!\n")

  # Summarize DE results
  de_results <- rctd_cside@de_results$res_gene_list
  if (length(de_results) > 0) {
    for (ct_name in names(de_results)) {
      sig_genes <- de_results[[ct_name]]
      n_sig <- nrow(sig_genes)
      cat(sprintf("    %-30s %d significant genes\n", ct_name, n_sig))
      if (n_sig > 0) {
        cat("      Top genes:\n")
        print(head(sig_genes, 5))
      }
    }
  }

  saveRDS(rctd_cside, out_cside)
  cat("  Saved:", out_cside, "\n")

}, error = function(e) {
  cat("  C-SIDE error:", conditionMessage(e), "\n")
  # Save error info so we can diagnose without re-running RCTD
  saveRDS(
    list(sample_id = sid, status = "cside_error", error = conditionMessage(e)),
    out_cside
  )
})

cat("\n=== Sample", sid, "done ===\n")
cat("Finished:", format(Sys.time()), "\n")
