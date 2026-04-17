#!/usr/bin/env Rscript
# ============================================================
# 05_merge_and_compare.R
#
# Run AFTER all array jobs complete. Loads individual sample
# RDS files, runs multi-sample comparisons, and generates
# summary statistics and plots.
# ============================================================

.libPaths(c(
  "/dartfs-hpc/rc/home/k/f006fpk/R/x86_64-redhat-linux-gnu-library/4.5",
  "/dartfs-hpc/admin/opt/el8/Rlibs/4.5",
  "/dartfs-hpc/admin/opt/el8/R/4.5.1/lib/R/library"
))

suppressPackageStartupMessages({
  library(SpatialCellChat)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ============================================================
# PATHS
# ============================================================
base_dir   <- "/dartfs/rc/nosnapshots/V/VaickusL-nb/EDIT_Students/users/alos/colon_work/adata"
meta_path  <- file.path(base_dir, "data_all.csv")
output_dir <- file.path(base_dir, "spatialcellchat_v3_results")
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# CLINICAL METADATA
# ============================================================
clinical <- read.csv(meta_path, stringsAsFactors = FALSE)
colnames(clinical) <- trimws(colnames(clinical))
clinical$deident <- as.character(clinical$deident)
clinical$MMR_status <- ifelse(clinical$MLH1 == 0 | clinical$PMS2 == 0, "dMMR", "pMMR")
clinical$met_group <- case_when(
  clinical$Distant_Mets == TRUE | clinical$Distant_Mets == "TRUE" ~ "Distant_Mets",
  clinical$ln_only == TRUE | clinical$ln_only == "TRUE" ~ "LN_Mets",
  TRUE ~ "No_Mets"
)

# ============================================================
# LOAD ALL COMPLETED SAMPLES
# ============================================================
rds_files <- list.files(output_dir, pattern = "^spatialcellchat_.*\\.rds$", full.names = TRUE)
rds_files <- rds_files[!grepl("all_samples", rds_files)]

cellchat_list <- list()
for (f in rds_files) {
  sid <- gsub("^spatialcellchat_|\\.rds$", "", basename(f))
  cellchat_list[[sid]] <- readRDS(f)
}

cat("Loaded", length(cellchat_list), "samples\n")

# Check for missing
prep_files <- list.files(file.path(base_dir, "prepped_samples"), pattern = "^prep_.*\\.rds$")
expected <- gsub("^prep_|\\.rds$", "", prep_files)
missing <- setdiff(expected, names(cellchat_list))
if (length(missing) > 0) {
  cat("WARNING: Missing samples:", paste(missing, collapse = ", "), "\n")
  writeLines(missing, file.path(output_dir, "missing_samples.txt"))
}

saveRDS(cellchat_list, file.path(output_dir, "spatialcellchat_all_samples.rds"))

# ============================================================
# MULTI-SAMPLE COMPARISON
# ============================================================
processed_ids <- names(cellchat_list)
clin_sub <- clinical[clinical$deident %in% processed_ids, , drop = FALSE]

compare_groups <- function(cc_list, a_ids, b_ids, group_names = c("A", "B")) {
  a_ids <- intersect(a_ids, names(cc_list))
  b_ids <- intersect(b_ids, names(cc_list))
  cat("Group", group_names[1], ":", length(a_ids), "samples\n")
  cat("Group", group_names[2], ":", length(b_ids), "samples\n")
  if (length(a_ids) < 2 || length(b_ids) < 2) stop("Need >= 2 per group")

  all_types <- unique(unlist(lapply(cc_list[c(a_ids, b_ids)],
                                    function(cc) levels(cc@idents))))
  la <- lapply(a_ids, function(s) liftCellChat(cc_list[[s]], group.new = all_types))
  lb <- lapply(b_ids, function(s) liftCellChat(cc_list[[s]], group.new = all_types))
  names(la) <- a_ids; names(lb) <- b_ids

  cc_a <- mergeCellChat(la, add.names = a_ids, cell.prefix = TRUE)
  cc_b <- mergeCellChat(lb, add.names = b_ids, cell.prefix = TRUE)
  mergeCellChat(list(cc_a, cc_b), add.names = group_names, cell.prefix = TRUE)
}

# dMMR vs pMMR
dmmr_ids <- clin_sub$deident[clin_sub$MMR_status == "dMMR"]
pmmr_ids <- clin_sub$deident[clin_sub$MMR_status == "pMMR"]

if (length(intersect(dmmr_ids, processed_ids)) >= 2 &&
    length(intersect(pmmr_ids, processed_ids)) >= 2) {
  tryCatch({
    cc_mmr <- compare_groups(cellchat_list, dmmr_ids, pmmr_ids, c("dMMR", "pMMR"))
    saveRDS(cc_mmr, file.path(output_dir, "comparison_dMMR_vs_pMMR.rds"))

    pdf(file.path(output_dir, "plots", "comparison_dMMR_vs_pMMR.pdf"), width = 12, height = 10)
    compareInteractions(cc_mmr, show.legend = TRUE, group = c(1, 2))
    compareInteractions(cc_mmr, show.legend = TRUE, group = c(1, 2), measure = "weight")
    dev.off()
    cat("✓ dMMR vs pMMR comparison complete\n")
  }, error = function(e) cat("✗ dMMR vs pMMR FAILED:", conditionMessage(e), "\n"))
}

# Distant mets vs no mets
distant_ids <- clin_sub$deident[clin_sub$met_group == "Distant_Mets"]
nomets_ids  <- clin_sub$deident[clin_sub$met_group == "No_Mets"]

if (length(intersect(distant_ids, processed_ids)) >= 2 &&
    length(intersect(nomets_ids, processed_ids)) >= 2) {
  tryCatch({
    cc_met <- compare_groups(cellchat_list, distant_ids, nomets_ids, c("Distant_Mets", "No_Mets"))
    saveRDS(cc_met, file.path(output_dir, "comparison_DistantMets_vs_NoMets.rds"))
    cat("✓ Distant mets vs no mets comparison complete\n")
  }, error = function(e) cat("✗ Mets comparison FAILED:", conditionMessage(e), "\n"))
}

# Right vs left/rectum
clin_sub$site_group <- case_when(
  grepl("Right|Cecum|Hepatic|Transverse", clin_sub$site, ignore.case = TRUE) ~ "Right_colon",
  grepl("Left|Sigmoid|Splenic", clin_sub$site, ignore.case = TRUE) ~ "Left_colon",
  grepl("Rectum", clin_sub$site, ignore.case = TRUE) ~ "Rectum",
  TRUE ~ "Other"
)
right_ids <- clin_sub$deident[clin_sub$site_group == "Right_colon"]
left_ids  <- clin_sub$deident[clin_sub$site_group %in% c("Left_colon", "Rectum")]

if (length(intersect(right_ids, processed_ids)) >= 2 &&
    length(intersect(left_ids, processed_ids)) >= 2) {
  tryCatch({
    cc_site <- compare_groups(cellchat_list, right_ids, left_ids, c("Right_colon", "Left_colon_Rectum"))
    saveRDS(cc_site, file.path(output_dir, "comparison_Right_vs_Left.rds"))
    cat("✓ Right vs left comparison complete\n")
  }, error = function(e) cat("✗ Site comparison FAILED:", conditionMessage(e), "\n"))
}

# ============================================================
# QUANTITATIVE SUMMARY
# ============================================================
summary_list <- lapply(names(cellchat_list), function(sid) {
  cc <- cellchat_list[[sid]]
  nw <- cc@net$weight; nc <- cc@net$count
  ex <- function(s, t, m) {
    if (s %in% rownames(m) && t %in% colnames(m)) m[s, t] else NA_real_
  }
  data.frame(
    sample_id = sid,
    n_spots = ncol(cc@data),
    n_cell_types = length(levels(cc@idents)),
    cell_types_present = paste(levels(cc@idents), collapse = ";"),
    total_interactions = sum(nc), total_strength = sum(nw),
    tumor_to_cd8 = ex("Tumor", "CD8_T", nw),
    tumor_to_cd8exh = ex("Tumor", "CD8_T_exhausted", nw),
    tumor_to_treg = ex("Tumor", "Treg", nw),
    tumor_to_cd4 = ex("Tumor", "CD4_T", nw),
    tumor_to_macrophage = ex("Tumor", "Macrophage", nw),
    tumor_to_dc = ex("Tumor", "DC", nw),
    caf_to_cd8 = ex("CAF", "CD8_T", nw),
    caf_to_cd8exh = ex("CAF", "CD8_T_exhausted", nw),
    caf_to_treg = ex("CAF", "Treg", nw),
    caf_to_macrophage = ex("CAF", "Macrophage", nw),
    caf_to_tumor = ex("CAF", "Tumor", nw),
    fibroblast_to_cd8 = ex("Fibroblast", "CD8_T", nw),
    myofibro_to_cd8 = ex("Myofibroblast", "CD8_T", nw),
    cd8_to_tumor = ex("CD8_T", "Tumor", nw),
    cd8exh_to_tumor = ex("CD8_T_exhausted", "Tumor", nw),
    macrophage_to_cd8 = ex("Macrophage", "CD8_T", nw),
    macrophage_to_treg = ex("Macrophage", "Treg", nw),
    dc_to_cd8 = ex("DC", "CD8_T", nw),
    treg_to_cd8 = ex("Treg", "CD8_T", nw),
    tumor_to_caf = ex("Tumor", "CAF", nw),
    tumor_to_endothelial = ex("Tumor", "Endothelial", nw),
    stringsAsFactors = FALSE
  )
})

interaction_summary <- do.call(rbind, summary_list)

interaction_clinical <- merge(
  interaction_summary,
  clinical[, c("deident", "Grade", "MMR_status", "met_group",
               "node_stage", "age", "sex", "site",
               "Distant_Mets", "Death_from_cancer")],
  by.x = "sample_id", by.y = "deident", all.x = TRUE
)

write.csv(interaction_clinical,
          file.path(output_dir, "interaction_summary_with_clinical.csv"),
          row.names = FALSE)

# ============================================================
# INTERACTION BOXPLOTS
# ============================================================
plot_cols <- c(
  "tumor_to_cd8", "tumor_to_cd8exh", "tumor_to_treg",
  "caf_to_cd8", "caf_to_tumor", "cd8_to_tumor",
  "caf_to_macrophage", "macrophage_to_cd8"
)

plot_df <- interaction_clinical %>%
  select(sample_id, MMR_status, all_of(plot_cols)) %>%
  pivot_longer(cols = all_of(plot_cols),
               names_to = "interaction", values_to = "strength") %>%
  filter(!is.na(strength), !is.na(MMR_status))

p1 <- ggplot(plot_df, aes(x = interaction, y = strength, fill = MMR_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.75)) +
  geom_jitter(aes(color = MMR_status),
              position = position_jitterdodge(jitter.width = 0.15),
              size = 1, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(title = "SpatialCellChat v3: Immune-Tumor-Fibroblast interactions (dMMR vs pMMR)",
       subtitle = paste0("N = ", nrow(interaction_summary), " samples | multi-label 0.05 threshold"),
       y = "Interaction strength", x = "") +
  scale_fill_manual(values = c("dMMR" = "#E41A1C", "pMMR" = "#377EB8")) +
  scale_color_manual(values = c("dMMR" = "#E41A1C", "pMMR" = "#377EB8"))

ggsave(file.path(output_dir, "plots", "interaction_boxplots_MMR.pdf"),
       plot = p1, width = 12, height = 7)

plot_df_met <- interaction_clinical %>%
  select(sample_id, met_group, all_of(plot_cols)) %>%
  pivot_longer(cols = all_of(plot_cols),
               names_to = "interaction", values_to = "strength") %>%
  filter(!is.na(strength), !is.na(met_group))

p2 <- ggplot(plot_df_met, aes(x = interaction, y = strength, fill = met_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.75)) +
  geom_jitter(aes(color = met_group),
              position = position_jitterdodge(jitter.width = 0.15),
              size = 1, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(title = "SpatialCellChat v3: Interactions by Metastatic Status",
       y = "Interaction strength", x = "") +
  scale_fill_manual(values = c("No_Mets" = "#00AFBB", "LN_Mets" = "#E7B800",
                               "Distant_Mets" = "#FC4E07")) +
  scale_color_manual(values = c("No_Mets" = "#00AFBB", "LN_Mets" = "#E7B800",
                                "Distant_Mets" = "#FC4E07"))

ggsave(file.path(output_dir, "plots", "interaction_boxplots_metastasis.pdf"),
       plot = p2, width = 12, height = 7)

# ============================================================
# STATISTICAL TESTS
# ============================================================
cat("\n=== Wilcoxon tests: dMMR vs pMMR ===\n")
for (col in plot_cols) {
  vals_dmmr <- interaction_clinical[[col]][interaction_clinical$MMR_status == "dMMR"]
  vals_pmmr <- interaction_clinical[[col]][interaction_clinical$MMR_status == "pMMR"]
  vals_dmmr <- vals_dmmr[!is.na(vals_dmmr)]
  vals_pmmr <- vals_pmmr[!is.na(vals_pmmr)]
  if (length(vals_dmmr) >= 2 && length(vals_pmmr) >= 2) {
    wt <- wilcox.test(vals_dmmr, vals_pmmr)
    cat(sprintf("  %-25s p = %.4f  (dMMR med: %.4f, pMMR med: %.4f)\n",
                col, wt$p.value, median(vals_dmmr), median(vals_pmmr)))
  } else {
    cat(sprintf("  %-25s SKIPPED (insufficient data)\n", col))
  }
}

cat("\n=== Kruskal-Wallis: Metastatic groups ===\n")
for (col in plot_cols) {
  groups <- split(interaction_clinical[[col]], interaction_clinical$met_group)
  groups <- lapply(groups, function(x) x[!is.na(x)])
  groups <- groups[sapply(groups, length) >= 2]
  if (length(groups) >= 2) {
    kt <- kruskal.test(groups)
    cat(sprintf("  %-25s p = %.4f\n", col, kt$p.value))
  }
}

# ============================================================
# PATHWAY SUMMARY
# ============================================================
pathway_summary <- lapply(names(cellchat_list), function(sid) {
  cc <- cellchat_list[[sid]]
  if (is.null(cc@netP$prob)) return(NULL)
  pathways <- dimnames(cc@netP$prob)[[3]]
  pw_data <- lapply(pathways, function(pw) {
    tryCatch({
      prob_mat <- cc@netP$prob[, , pw]
      data.frame(sample_id = sid, pathway = pw,
                 total_prob = sum(prob_mat, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
  })
  do.call(rbind, pw_data)
})

pathway_df <- do.call(rbind, pathway_summary)

top_pathways <- pathway_df %>%
  group_by(pathway) %>%
  summarise(
    mean_prob = mean(total_prob, na.rm = TRUE),
    median_prob = median(total_prob, na.rm = TRUE),
    n_samples_active = sum(total_prob > 0),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_prob)) %>%
  head(20)

cat("\nTop 20 active pathways:\n")
print(top_pathways)

write.csv(pathway_df, file.path(output_dir, "pathway_summary_all_samples.csv"),
          row.names = FALSE)

# ============================================================
# SAVE
# ============================================================
saveRDS(interaction_summary, file.path(output_dir, "interaction_summary.rds"))
saveRDS(interaction_clinical, file.path(output_dir, "interaction_clinical.rds"))
saveRDS(pathway_df, file.path(output_dir, "pathway_summary.rds"))

cat("\nAll results saved to:", output_dir, "\n")
cat("Files:\n")
print(list.files(output_dir, recursive = TRUE))
