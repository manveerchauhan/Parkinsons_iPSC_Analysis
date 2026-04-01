#!/usr/bin/env Rscript
# =============================================================================
# sample_makeup.R
# =============================================================================
# Visualise the donor × condition composition of the filtered PPMI iPSC
# scRNA-seq dataset after filter_and_annotate.R.
#
# Adapted from: ch3_cellTypeProportion_script.R (Neurodevelopmental_Models_analysis)
# Author : Manveer Chauhan
# Usage  : Rscript sample_makeup.R  (run from Parkinsons_iPSC_Analysis/)
#
# Input  : export/seurat_filtered.rds
# Output : output_files/cell_type_proportions/iPSC_condition_sampleMakeup.pdf
#          output_files/cell_type_proportions/iPSC_donor_sampleMakeup.pdf
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_RDS  <- "export/seurat_qc_filtered.rds"
OUTPUT_DIR <- "output_files/cell_type_proportions"

# =============================================================================
# LIBRARIES
# =============================================================================
library(Seurat)
library(tidyverse)
library(gridExtra)

# Use ggmin theme if available, otherwise fall back to theme_bw
if (requireNamespace("ggmin", quietly = TRUE)) {
  theme_set(ggmin::theme_min())
} else {
  theme_set(theme_bw())
}

# =============================================================================
# HELPER
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
# EXTRACT METADATA
# =============================================================================
section("EXTRACTING METADATA")

DONOR_COL <- "PPMI Donor ID"

meta <- seurat_obj@meta.data %>%
  select(donor = all_of(DONOR_COL), condition) %>%
  group_by(donor, condition) %>%
  summarise(n_cells = n(), .groups = "drop")

cat("\n  Cell counts per donor × condition:\n")
print(as.data.frame(meta))

cat("\n  Cell counts per condition:\n")
print(meta %>% group_by(condition) %>% summarise(total_cells = sum(n_cells)))

# =============================================================================
# PLOT 1 — CONDITION-LEVEL (global)
# Two panels: stacked counts + stacked proportions
# x = condition, fill = donor
# =============================================================================
section("PLOT 1: CONDITION-LEVEL COMPOSITION")

# Colour palette — one colour per donor (10 donors)
n_donors <- length(unique(meta$donor))
donor_colours <- setNames(
  scales::hue_pal()(n_donors),
  sort(unique(meta$donor))
)

# Count version
p_cond_count <- ggplot(meta, aes(x = condition, y = n_cells, fill = donor)) +
  geom_bar(stat = "identity", position = "stack", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = donor_colours, name = "PPMI Donor ID") +
  labs(title = "Cells per condition (stacked by donor)",
       x = "Condition", y = "Cell count") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Proportion version
p_cond_prop <- ggplot(meta, aes(x = condition, y = n_cells, fill = donor)) +
  geom_bar(stat = "identity", position = "fill", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = donor_colours, name = "PPMI Donor ID") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Donor proportions per condition",
       x = "Condition", y = "Proportion of cells") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

cond_grid <- grid.arrange(p_cond_count, p_cond_prop, ncol = 2)

# =============================================================================
# PLOT 2 — DONOR-LEVEL (replicate)
# x = donor, y = cell count, fill = condition, faceted by condition
# =============================================================================
section("PLOT 2: DONOR-LEVEL COMPOSITION")

condition_colours <- c(
  "Idiopathic_PD" = "#E64B35",
  "Control"       = "#4DBBD5"
)

p_donor <- ggplot(meta, aes(x = donor, y = n_cells, fill = condition)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = condition_colours, name = "Condition") +
  facet_wrap(~ condition, scales = "free_x") +
  labs(title = "Cells per donor (faceted by condition)",
       x = "PPMI Donor ID", y = "Cell count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

# =============================================================================
# SAVE
# =============================================================================
section("SAVING PLOTS")

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("  Created output directory:", OUTPUT_DIR, "\n")
}

out_cond  <- file.path(OUTPUT_DIR, "iPSC_condition_sampleMakeup.pdf")
out_donor <- file.path(OUTPUT_DIR, "iPSC_donor_sampleMakeup.pdf")

ggsave(out_cond,  cond_grid, width = 12, height = 6, dpi = 300)
ggsave(out_donor, p_donor,   width = 14, height = 6, dpi = 300)

cat("  Saved:", out_cond,  "\n")
cat("  Saved:", out_donor, "\n")

cat(sprintf("
  INPUT  : %s
  OUTPUT : %s
           %s
  Total cells plotted: %d

", INPUT_RDS, out_cond, out_donor, sum(meta$n_cells)))
