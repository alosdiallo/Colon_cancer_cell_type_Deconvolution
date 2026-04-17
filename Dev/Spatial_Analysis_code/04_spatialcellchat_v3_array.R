#!/usr/bin/env Rscript
# ============================================================
# 04_spatialcellchat_v3_array.R
#
# Processes ONE sample based on the SLURM array task ID.
# Called by submit_sccv3_array.sh
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 04_spatialcellchat_v3_array.R <TASK_ID>")
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
  library(SpatialCellChat)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

# ============================================================
# PARALLEL SETTINGS
# ============================================================
nworkers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(nworkers) || nworkers < 1) nworkers <- 1

setEnvironment(
  workers = nworkers,
  future.globals.maxSize = 100000 * 1024^2,
  conda_env = NULL
)
cat("Workers:", nworkers, "\n")

# ============================================================
# PATHS
# ============================================================
base_dir   <- "/dartfs/rc/nosnapshots/V/VaickusL-nb/EDIT_Students/users/alos/colon_work/adata"
prep_dir   <- file.path(base_dir, "prepped_samples")
output_dir <- file.path(base_dir, "spatialcellchat_v3_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# FIND SAMPLE FOR THIS TASK
# ============================================================
prep_files <- sort(list.files(prep_dir, pattern = "^prep_.*\\.rds$", full.names = FALSE))
sample_ids <- gsub("^prep_|\\.rds$", "", prep_files)

if (task_id >= length(sample_ids)) {
  cat("Task ID", task_id, "exceeds number of samples (", length(sample_ids), "). Exiting.\n")
  quit(status = 0)
}

sid <- sample_ids[task_id + 1]  # +1 because R is 1-indexed, array is 0-indexed
cat("Task ID:", task_id, "| Sample:", sid, "\n")

# Skip if already completed
out_path <- file.path(output_dir, paste0("spatialcellchat_", sid, ".rds"))
if (file.exists(out_path)) {
  cat("Sample", sid, "already complete. Skipping.\n")
  quit(status = 0)
}

# ============================================================
# CELL TYPE MAPPING (needed for rebuild)
# ============================================================
map_c2l_to_celltype <- function(col_name) {
  clean <- gsub("meanscell_abundance_w_sf_", "", col_name)
  clean <- gsub("means_cell_abundance_w_sf_", "", clean)
  if (grepl("^Tumor", clean)) return("Tumor")
  if (grepl("^cE", clean)) return("Epithelial")
  if (grepl("^cS", clean)) {
    num <- as.integer(gsub("\\D", "", sub("\\s.*", "", clean)))
    if (!is.na(num)) {
      if (num >= 1 & num <= 14) return("Endothelial")
      if (num == 32) return("Smooth_muscle")
      if (num %in% c(27:31)) return("CAF")
      if (num == 26) return("Myofibroblast")
      if (num %in% c(22, 25)) return("Fibroblast")
      if (num %in% c(15:20)) return("Pericyte")
      if (num == 33) return("Schwann")
    }
    return("Fibroblast")
  }
  if (grepl("^cB", clean)) return("B_cell")
  if (grepl("^cP", clean)) return("Plasma")
  if (grepl("^cMA", clean)) return("Mast")
  if (grepl("^cM", clean)) {
    if (grepl("cM10", clean)) return("Granulocyte")
    if (grepl("cM01|cM02", clean)) return("Macrophage")
    return("DC")
  }
  if (grepl("^cTNI", clean)) {
    num <- as.integer(gsub("\\D", "", sub("\\s.*", "", clean)))
    if (!is.na(num)) {
      if (num >= 1 & num <= 7) return("CD4_T")
      if (num %in% c(8, 9)) return("Treg")
      if (num %in% c(10:13)) return("CD8_T")
      if (num %in% c(14:16)) return("CD8_T_exhausted")
      if (num %in% c(17, 18)) return("gdT")
      if (num %in% c(19:22)) return("Other_T")
      if (num %in% c(23:25)) return("NK")
      if (num == 26) return("ILC")
    }
    return("Other_T")
  }
  return("Other")
}

# ============================================================
# MULTI-LABEL REBUILD
# ============================================================
rebuild_multilabel <- function(prep, prop_threshold = 0.05,
                                max_spots_per_type = 500) {

  cat("  Rebuilding with multi-label (threshold:", prop_threshold, ")\n")

  data_norm <- tryCatch(
    GetAssayData(prep$seu, layer = "data", assay = "RNA"),
    error = function(e) GetAssayData(prep$seu, slot = "data", assay = "RNA")
  )
  counts_raw <- tryCatch(
    GetAssayData(prep$seu, layer = "counts", assay = "RNA"),
    error = function(e) GetAssayData(prep$seu, slot = "counts", assay = "RNA")
  )

  c2l_agg <- prep$c2l_agg
  row_totals <- rowSums(c2l_agg)
  row_totals[row_totals == 0] <- 1
  c2l_prop <- c2l_agg / row_totals

  orig_coords <- as.data.frame(prep$coords)
  orig_meta   <- prep$meta

  new_labels  <- c()
  new_spots   <- c()
  spot_origins <- c()

  for (i in seq_len(nrow(c2l_prop))) {
    spot_name <- rownames(c2l_prop)[i]
    above <- names(which(c2l_prop[i, ] >= prop_threshold))
    dominant <- names(which.max(c2l_agg[i, ]))
    above <- unique(c(above, dominant))
    above <- above[!above %in% c("Other", "Schwann", "ILC", "Other_T")]

    for (ct in above) {
      new_id <- paste0(spot_name, "_", ct)
      new_labels  <- c(new_labels, ct)
      new_spots   <- c(new_spots, new_id)
      spot_origins <- c(spot_origins, spot_name)
    }
  }

  cat("  Original spots:", nrow(c2l_prop), "\n")
  cat("  Multi-label spots:", length(new_spots), "\n")
  cat("  Label distribution:\n")
  print(sort(table(new_labels), decreasing = TRUE))

  set.seed(42)
  keep_idx <- c()
  for (ct in unique(new_labels)) {
    ct_idx <- which(new_labels == ct)
    if (length(ct_idx) > max_spots_per_type) {
      ct_idx <- sample(ct_idx, max_spots_per_type)
    }
    keep_idx <- c(keep_idx, ct_idx)
  }
  keep_idx <- sort(keep_idx)

  new_labels   <- new_labels[keep_idx]
  new_spots    <- new_spots[keep_idx]
  spot_origins <- spot_origins[keep_idx]

  cat("  After downsampling:", length(new_spots), "spots\n")
  cat("  Post-downsample distribution:\n")
  print(sort(table(new_labels), decreasing = TRUE))

  dup_counts <- counts_raw[, spot_origins, drop = FALSE]
  colnames(dup_counts) <- new_spots

  dup_coords <- orig_coords[spot_origins, , drop = FALSE]
  rownames(dup_coords) <- new_spots

  dup_meta <- data.frame(
    labels = new_labels,
    histology = orig_meta$histology[match(spot_origins, rownames(orig_meta))],
    row.names = new_spots
  )

  dup_c2l_prop <- c2l_prop[spot_origins, , drop = FALSE]
  rownames(dup_c2l_prop) <- new_spots
  rs <- rowSums(dup_c2l_prop); rs[rs == 0] <- 1
  dup_c2l_prop <- dup_c2l_prop / rs

  seu <- CreateSeuratObject(counts = dup_counts, meta.data = dup_meta)
  seu <- NormalizeData(seu)

  list(
    seu             = seu,
    coords          = dup_coords,
    spatial_factors = list(ratio = 0.5686331, tol = 65 / 2),
    sample_id       = prep$sample_id,
    c2l_proportions = dup_c2l_prop
  )
}

# ============================================================
# RUN SPATIALCELLCHAT
# ============================================================
run_spatialcellchat <- function(prep, min_cells = 10) {

  cat("\n--- SpatialCellChat for:", prep$sample_id, "---\n")

  data.input <- tryCatch(
    GetAssayData(prep$seu, layer = "data", assay = "RNA"),
    error = function(e) GetAssayData(prep$seu, slot = "data", assay = "RNA")
  )
  meta.input <- prep$seu@meta.data

  chat <- createSpatialCellChat(
    object = data.input,
    meta = meta.input,
    group.by = "labels",
    datatype = "spatial",
    coordinates = prep$coords,
    spatial.factors = prep$spatial_factors
  )

  CellChatDB.use <- subsetDB(
    CellChatDB.human,
    search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
    non_protein = FALSE
  )
  chat@DB <- CellChatDB.use

  chat <- subsetData(chat)
  chat <- preProcessing(chat)
  chat <- identifyOverExpressedGenes(chat, selection.method = "meringue", do.grid = FALSE)
  chat <- identifyOverExpressedInteractions(chat, variable.both = FALSE)

  chat <- computeCommunProb(
    chat, distance.use = TRUE, scale.distance = 0.2,
    contact.dependent = TRUE,
    interaction.range = 250, contact.range = 100
  )

  chat <- filterProbability(chat)
  chat <- filterCommunication(chat, min.cells = NULL, min.links = 10, min.cells.sr = 10)

  # Group-level using C2L decomposition
  chat_groups <- levels(chat@idents)
  c2l <- prep$c2l_proportions
  shared <- intersect(colnames(c2l), chat_groups)
  c2l <- c2l[, shared, drop = FALSE]
  missing <- setdiff(chat_groups, shared)
  if (length(missing) > 0) {
    zero_mat <- matrix(0, nrow = nrow(c2l), ncol = length(missing),
                       dimnames = list(rownames(c2l), missing))
    c2l <- cbind(c2l, zero_mat)
  }
  c2l <- c2l[, chat_groups, drop = FALSE]
  rs <- rowSums(c2l); rs[rs == 0] <- 1; c2l <- c2l / rs

  chat <- computeAvgCommunProb_Visium(
    chat, cell.type.decomposition = c2l,
    avg.type = "avg", nboot = 50, do.permutation = TRUE
  )
  chat <- filterCommunication(chat, min.cells = min_cells, min.links = NULL, min.cells.sr = NULL)

  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)

  chat <- netAnalysis_computeCentrality(chat, slot.name = "net", do.group = FALSE, degree.only = TRUE)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "net", do.group = TRUE, degree.only = TRUE)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP", do.group = FALSE, degree.only = TRUE)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP", do.group = TRUE, degree.only = TRUE)

  cat("  Cell types:", paste(levels(chat@idents), collapse = ", "), "\n")
  cat("  Active pathways:", length(chat@netP$pathways), "\n")

  return(chat)
}

# ============================================================
# PROCESS THIS SAMPLE
# ============================================================
tryCatch({
  cat("\n========================================\n")
  cat("Processing sample:", sid, "\n")
  cat("========================================\n")

  prep <- readRDS(file.path(prep_dir, paste0("prep_", sid, ".rds")))
  prep_ml <- rebuild_multilabel(prep, prop_threshold = 0.05, max_spots_per_type = 500)
  cc <- run_spatialcellchat(prep_ml, min_cells = 10)

  saveRDS(cc, out_path)
  cat("✓ Sample", sid, "complete\n")

}, error = function(e) {
  cat("✗ Sample", sid, "FAILED:", conditionMessage(e), "\n")
  quit(status = 1)
})

cat("Done:", format(Sys.time()), "\n")
