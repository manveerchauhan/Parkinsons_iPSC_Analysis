#!/usr/bin/env Rscript
# =============================================================================
# qc_visualisation.R
# =============================================================================
# Per-donor QC metric visualisation for PPMI iPSC scRNA-seq data.
# Purpose: verify dataset quality BEFORE Harmony integration (visualisation only
# — no filtering; PPMI cells were already QC'd in the upstream published pipeline).
#
# Author : Manveer Chauhan
# Usage  : Rscript qc_visualisation.R  (run from Parkinsons_iPSC_Analysis/)
#
# Input  : export/seurat_filtered.rds
# Output : output_files/qc/iPSC_QC_per_donor.pdf
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_RDS  <- "export/seurat_qc_filtered.rds"
OUTPUT_DIR <- "output_files/qc"
OUTPUT_PDF <- file.path(OUTPUT_DIR, "iPSC_QC_per_donor.pdf")
DONOR_COL  <- "PPMI Donor ID"

# =============================================================================
# LIBRARIES
# =============================================================================
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)

if (requireNamespace("ggmin", quietly = TRUE)) {
  theme_set(ggmin::theme_min())
} else {
  theme_set(theme_bw())
}

# Condition colour palette (consistent with rest of pipeline)
condition_colours <- c("Idiopathic_PD" = "#E64B35", "Control" = "#4DBBD5")

# =============================================================================
# HELPERS
# =============================================================================
section <- function(title) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  ", title, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

# =============================================================================
# LOAD
# =============================================================================
section("LOADING SEURAT OBJECT")

if (!file.exists(INPUT_RDS)) {
  stop("Could not find: ", INPUT_RDS, "\n  Run filter_and_annotate.R first.")
}
seurat_obj <- readRDS(INPUT_RDS)
cat("  Loaded:", INPUT_RDS, "\n")
cat("  Total cells:", ncol(seurat_obj), "\n")

# =============================================================================
# CALCULATE PERCENT MITOCHONDRIA
# =============================================================================
section("CALCULATING PERCENT MITOCHONDRIA")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
cat("  percent.mt added (genes matching ^MT-)\n")
cat("  Median percent.mt:", round(median(seurat_obj$percent.mt), 2), "%\n")

# =============================================================================
# SUMMARY TABLE
# =============================================================================
section("PER-DONOR QC SUMMARY")

meta <- seurat_obj@meta.data %>%
  rename(donor = all_of(DONOR_COL))

summary_tbl <- meta %>%
  group_by(donor, condition) %>%
  summarise(
    n_cells       = n(),
    mean_nFeature = round(mean(nFeature_RNA), 0),
    med_nFeature  = round(median(nFeature_RNA), 0),
    mean_nCount   = round(mean(nCount_RNA), 0),
    med_nCount    = round(median(nCount_RNA), 0),
    mean_pct_mt   = round(mean(percent.mt), 2),
    med_pct_mt    = round(median(percent.mt), 2),
    .groups       = "drop"
  ) %>%
  arrange(condition, donor)

cat("\n")
print(as.data.frame(summary_tbl))

# =============================================================================
# PLOT HELPERS
# =============================================================================

# Ordered donor factor (group by condition for visual clarity)
donor_order <- summary_tbl %>%
  arrange(condition, donor) %>%
  pull(donor)

meta$donor    <- factor(meta$donor, levels = donor_order)
meta$condition <- factor(meta$condition, levels = c("Control", "Idiopathic_PD"))

# =============================================================================
# PAGE 1: VIOLIN PLOTS — nFeature, nCount, percent.mt per donor
# =============================================================================
p_nfeature <- ggplot(meta, aes(x = donor, y = nFeature_RNA, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3, alpha = 0.7) +
  scale_fill_manual(values = condition_colours) +
  labs(title = "Genes per cell (nFeature_RNA)",
       x = NULL, y = "nFeature_RNA", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ncount <- ggplot(meta, aes(x = donor, y = nCount_RNA, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3, alpha = 0.7) +
  scale_fill_manual(values = condition_colours) +
  labs(title = "UMI counts per cell (nCount_RNA)",
       x = NULL, y = "nCount_RNA", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_pctmt <- ggplot(meta, aes(x = donor, y = percent.mt, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3, alpha = 0.7) +
  scale_fill_manual(values = condition_colours) +
  labs(title = "Mitochondrial gene % (percent.mt)",
       x = "PPMI Donor ID", y = "percent.mt (%)", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

page1 <- grid.arrange(
  p_nfeature, p_ncount, p_pctmt,
  nrow = 3,
  top = textGrob("iPSC QC per Donor — nFeature / nCount / percent.mt",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# =============================================================================
# PAGE 2: SCATTER — nCount vs nFeature, faceted by donor, coloured by percent.mt
# =============================================================================

# Subsample per donor for readability (max 500 cells per donor in scatter)
set.seed(42)
meta_sub <- meta %>%
  group_by(donor) %>%
  slice_sample(prop = 1) %>%   # random shuffle within each donor
  slice_head(n = 500) %>%      # take up to 500 (all if donor has fewer)
  ungroup()

p_scatter_feat <- ggplot(meta_sub,
                          aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_viridis_c(name = "% mito", option = "magma") +
  facet_wrap(~ donor, ncol = 5) +
  labs(title = "nCount vs nFeature per donor (coloured by % mito)",
       subtitle = "Subsampled to 500 cells/donor for clarity",
       x = "nCount_RNA", y = "nFeature_RNA") +
  theme(strip.text = element_text(size = 7),
        axis.text  = element_text(size = 6))

# =============================================================================
# PAGE 3: SCATTER — nCount vs percent.mt, faceted by donor
# =============================================================================
p_scatter_mt <- ggplot(meta_sub,
                        aes(x = nCount_RNA, y = percent.mt, colour = condition)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_manual(values = condition_colours) +
  facet_wrap(~ donor, ncol = 5) +
  labs(title = "nCount vs percent.mt per donor",
       subtitle = "High percent.mt with low nCount → low-quality / dying cells",
       x = "nCount_RNA", y = "percent.mt (%)") +
  theme(strip.text = element_text(size = 7),
        axis.text  = element_text(size = 6))

# =============================================================================
# PAGE 4: CELL COUNT BAR — cells per donor, coloured by condition
# =============================================================================
count_df <- summary_tbl %>%
  mutate(donor = factor(donor, levels = donor_order))

p_cellcount <- ggplot(count_df, aes(x = donor, y = n_cells, fill = condition)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  geom_text(aes(label = scales::comma(n_cells)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = condition_colours) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Cell counts per donor",
       x = "PPMI Donor ID", y = "Cell count", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =============================================================================
# SAVE PDF
# =============================================================================
section("SAVING PDF")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

pdf(OUTPUT_PDF, width = 16, height = 10)
grid.draw(page1)
print(p_scatter_feat)
print(p_scatter_mt)
print(p_cellcount)
dev.off()

cat("  Saved:", OUTPUT_PDF, "\n")

cat(sprintf("
  INPUT  : %s
  OUTPUT : %s
  Cells  : %d across %d donors

  Per-donor summary printed above.
  Review PDF for:
    - Donors with unusually low nFeature (empty droplets / low quality)
    - Donors with high percent.mt (stressed / dying cells)
    - Donors with outlier cell counts relative to peers

", INPUT_RDS, OUTPUT_PDF, ncol(seurat_obj), length(unique(meta$donor))))
