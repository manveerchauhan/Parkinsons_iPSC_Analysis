#!/usr/bin/env Rscript
# =============================================================================
# iPSC_harmony_integration.R
# =============================================================================
# Harmony integration pipeline for PPMI Parkinson's iPSC scRNA-seq data.
#
# Adapted from: ch3_org_integration_harmony.R (Neurodevelopmental_Models_analysis)
# Author : Manveer Chauhan
# Usage  : Rscript iPSC_harmony_integration.R  (run from Parkinsons_iPSC_Analysis/)
#
# Input  : export/seurat_filtered.rds  (10 donors: 5 Idiopathic_PD, 5 Control)
# Output : output_files/integrated_objects/iPSC_PD_integrated_harmony.rds
#          output_files/integration_QC/iPSC_PD_integration_report.pdf
#
# Key differences from the organoid script:
#   - Single pre-merged object (no per-timepoint merge loop)
#   - Harmony batch variable: "PPMI Donor ID" (donor-level correction)
#   - Cell cycle scoring added (not pre-run in upstream QC)
#   - PDF includes condition + Timepoint (DIV stage) UMAPs as diagnostics
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_RDS       <- "export/seurat_qc_filtered.rds"
OUTPUT_RDS      <- "output_files/integrated_objects/iPSC_PD_integrated_harmony.rds"
OUTPUT_PDF      <- "output_files/integration_QC/iPSC_PD_integration_report.pdf"
DONOR_COL       <- "PPMI Donor ID"
CUSTOM_RES      <- NA    # NA = use silhouette-optimal resolution

# =============================================================================
# LIBRARIES
# =============================================================================
library(Seurat)
library(tidyverse)
library(harmony)
library(clustree)
library(gridExtra)
library(grid)

# =============================================================================
# EXTERNAL DEPENDENCIES (all local — self-contained repo)
# =============================================================================
load("cycle.rda")       # provides s.genes and g2m.genes
source("silhouette.R")  # provides optimize_silhouette()
source("qc_functions.R")# provides quantitative_elbow()

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================
dir.create("output_files/integrated_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("output_files/integration_QC",     recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("HARMONY INTEGRATION PIPELINE — PPMI iPSC\n")
cat("========================================\n\n")

# =============================================================================
# STEP 1: Load filtered Seurat object
# =============================================================================
cat("=== STEP 1: Loading filtered Seurat object ===\n")

if (!file.exists(INPUT_RDS)) {
  stop("Could not find: ", INPUT_RDS, "\n  Run filter_and_annotate.R first.")
}
seurat_obj <- readRDS(INPUT_RDS)
cat("  Loaded:", INPUT_RDS, "\n")
cat("  Total cells:", ncol(seurat_obj), "\n")

cat("\n  Cells per condition:\n")
print(table(Condition = seurat_obj@meta.data$condition))
cat("\n  Cells per donor:\n")
print(table(Donor = seurat_obj@meta.data[[DONOR_COL]]))

# =============================================================================
# STEP 2: Normalization and cell cycle scoring (must run on single-layer object)
# CellCycleScoring calls GetAssayData internally — fails on multi-layer v5 assay
# =============================================================================
cat("\n=== STEP 2: Normalization and cell cycle scoring (pre-split) ===\n")
# Join any split layers from merge() before normalizing (CellCycleScoring requires single layer)
seurat_obj <- JoinLayers(seurat_obj)
cat("  Layers joined (pre-normalization)\n")
seurat_obj <- NormalizeData(seurat_obj)
cat("  Normalization complete\n")

seurat_obj <- CellCycleScoring(seurat_obj,
                                s.features   = s_genes,
                                g2m.features = g2m_genes,
                                set.ident    = FALSE)
cat("  Cell cycle scoring complete\n")
cat("  Phase distribution:\n")
print(table(seurat_obj$Phase))

# Seurat v5: split RNA layers by donor AFTER scoring (IntegrateLayers requirement)
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]],
                              f = seurat_obj@meta.data[[DONOR_COL]])
cat("  RNA layers split by donor (Seurat v5 IntegrateLayers requirement)\n")
cat("  Layers:", paste(Layers(seurat_obj), collapse = ", "), "\n")

# =============================================================================
# STEP 3: Pre-integration processing (HVG → Scale → PCA → elbow)
# Normalization and cell cycle scoring already done above (pre-split)
# =============================================================================
pre_integration_processing <- function(seurat_obj, label = "iPSC_PD") {

  cat("\n=== Pre-integration processing for", label, "===\n")

  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = "vst",
                                     nfeatures = 2000)
  cat("  Found 2000 variable features\n")

  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  cat("  Data scaling complete\n")

  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  cat("  PCA complete\n")

  # Quantitative elbow to determine optimal PCs (+5 adjustment)
  elbow_result <- quantitative_elbow(seurat_obj, label, max_pcs = 50, pc_adjustment = 5)
  optimal_pcs  <- elbow_result$optimal_pcs_used

  return(list(
    seurat_obj  = seurat_obj,
    optimal_pcs = optimal_pcs,
    elbow_plots = elbow_result$plots
  ))
}

pre_integration <- pre_integration_processing(seurat_obj)
cat("\n  Optimal PCs (with +5 adjustment):", pre_integration$optimal_pcs, "\n")

# =============================================================================
# STEP 4: Harmony integration
# =============================================================================
harmony_integration <- function(seurat_obj, label, optimal_pcs, donor_col) {

  cat("\n=== Running Harmony integration for", label, "===\n")

  # Pre-integration UMAP for comparison
  seurat_obj <- RunUMAP(seurat_obj,
                        dims = 1:optimal_pcs,
                        reduction = "pca",
                        reduction.name = "umap.unintegrated")
  cat("  Created pre-integration UMAP\n")

  # Harmony integration — correct for donor-level batch effects
  seurat_obj <- IntegrateLayers(
    object         = seurat_obj,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    group.by       = donor_col,
    verbose        = TRUE
  )
  cat("  Harmony integration complete\n")

  # Collapse split layers back into a single layer for downstream tools
  seurat_obj <- JoinLayers(seurat_obj)
  cat("  Layers joined post-integration\n")

  # Post-integration UMAP
  seurat_obj <- RunUMAP(seurat_obj,
                        reduction = "harmony",
                        dims = 1:optimal_pcs,
                        reduction.name = "umap.harmony")
  cat("  Created post-integration UMAP\n")

  before_plot <- DimPlot(seurat_obj, reduction = "umap.unintegrated",
                         group.by = donor_col) +
    ggtitle(paste0(label, ": Before Integration")) +
    labs(color = "Donor")

  after_plot <- DimPlot(seurat_obj, reduction = "umap.harmony",
                        group.by = donor_col) +
    ggtitle(paste0(label, ": After Harmony Integration")) +
    labs(color = "Donor")

  return(list(
    seurat_obj = seurat_obj,
    plots      = list(before = before_plot, after = after_plot)
  ))
}

harmony_result <- harmony_integration(pre_integration$seurat_obj,
                                      "iPSC_PD",
                                      pre_integration$optimal_pcs,
                                      DONOR_COL)

# =============================================================================
# STEP 5: Post-integration clustering
# =============================================================================
post_integration_clustering <- function(seurat_obj, label, optimal_pcs,
                                        custom_resolution = NA) {

  cat("\n=== Post-integration clustering for", label, "===\n")

  seurat_obj <- FindNeighbors(seurat_obj,
                               reduction = "harmony",
                               dims = 1:optimal_pcs)
  cat("  FindNeighbors complete (harmony, dims 1:", optimal_pcs, ")\n")

  # Silhouette uses stats::dist() which is O(n^2) memory — subsample to avoid OOM
  # 5K cells → ~200 MB distance matrix; representative of full cluster structure
  cat("  Running silhouette analysis (subsampled to 5000 cells)...\n")
  MAX_SIL_CELLS <- 5000
  if (ncol(seurat_obj) > MAX_SIL_CELLS) {
    set.seed(42)
    sil_cells <- sample(colnames(seurat_obj), MAX_SIL_CELLS)
    seurat_sub <- seurat_obj[, sil_cells]
  } else {
    seurat_sub <- seurat_obj
  }
  sil_results <- optimize_silhouette(sobject      = seurat_sub,
                                     test_res     = seq(0.1, 1.2, by = 0.1),
                                     summary_plot = TRUE,
                                     reduction    = "harmony")

  sil_df <- as.data.frame(sil_results) %>%
    dplyr::rename(avg_sil_vals = sil_vals) %>%
    group_by(num_clusters) %>%
    dplyr::slice(which.max(avg_sil_vals)) %>%
    ungroup() %>%
    arrange(desc(avg_sil_vals))

  optimal_res <- sil_df$res_vals[1]
  cat("  Calculated optimal resolution:", optimal_res, "(silhouette:",
      round(sil_df$avg_sil_vals[1], 4), ")\n")

  if (!is.na(custom_resolution)) {
    final_res         <- custom_resolution
    resolution_method <- "CUSTOM OVERRIDE"
    cat("  Using CUSTOM resolution:", final_res,
        "(calculated optimal was:", optimal_res, ")\n")
  } else {
    final_res         <- optimal_res
    resolution_method <- "calculated"
    cat("  Using calculated optimal resolution:", final_res, "\n")
  }

  # All resolutions for clustree
  seurat_obj <- FindClusters(seurat_obj, resolution = seq(0.1, 1.2, by = 0.1))

  clustree_plot <- clustree(seurat_obj) +
    ggtitle(paste0(label, ": Clustree")) +
    labs(subtitle = paste0("Final resolution: ", final_res,
                           " (", resolution_method, ")"))

  # Final clustering
  seurat_obj <- FindClusters(seurat_obj, resolution = final_res)

  cluster_umap <- DimPlot(seurat_obj, reduction = "umap.harmony", label = TRUE) +
    ggtitle(paste0(label, ": Final Clusters")) +
    labs(subtitle = paste0("Resolution: ", final_res, " (", resolution_method, ")"),
         color = "Cluster")

  return(list(
    seurat_obj        = seurat_obj,
    optimal_res       = optimal_res,
    final_res         = final_res,
    resolution_method = resolution_method,
    sil_results       = sil_df,
    plots             = list(clustree = clustree_plot, cluster_umap = cluster_umap)
  ))
}

clustering_result <- post_integration_clustering(harmony_result$seurat_obj,
                                                  "iPSC_PD",
                                                  pre_integration$optimal_pcs,
                                                  custom_resolution = CUSTOM_RES)

# =============================================================================
# STEP 6: QC plots (cell cycle, nFeature, nCount)
# =============================================================================
generate_qc_plots <- function(seurat_obj, label) {

  cat("\n=== Generating QC plots for", label, "===\n")

  cellcycle_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                            group.by = "Phase") +
    ggtitle(paste0(label, ": Cell Cycle Phase")) +
    labs(color = "Phase")

  nfeature_umap <- FeaturePlot(seurat_obj, reduction = "umap.harmony",
                               features = "nFeature_RNA") +
    ggtitle(paste0(label, ": nFeature_RNA"))

  ncount_umap <- FeaturePlot(seurat_obj, reduction = "umap.harmony",
                             features = "nCount_RNA") +
    ggtitle(paste0(label, ": nCount_RNA"))

  return(list(cellcycle = cellcycle_umap,
              nfeature  = nfeature_umap,
              ncount    = ncount_umap))
}

qc_plots <- generate_qc_plots(clustering_result$seurat_obj, "iPSC_PD")

# =============================================================================
# STEP 7: Extra diagnostic plots
# =============================================================================
cat("\n=== Generating extra diagnostic plots ===\n")

final_obj <- clustering_result$seurat_obj

donor_umap <- DimPlot(final_obj, reduction = "umap.harmony",
                      group.by = DONOR_COL) +
  ggtitle("iPSC_PD: PPMI Donor ID") +
  labs(color = "Donor")

condition_umap <- DimPlot(final_obj, reduction = "umap.harmony",
                          group.by = "condition") +
  ggtitle("iPSC_PD: Condition") +
  labs(color = "Condition")

before_condition_umap <- DimPlot(final_obj, reduction = "umap.unintegrated",
                                  group.by = "condition") +
  ggtitle("iPSC_PD: Before Integration (Condition)") +
  labs(color = "Condition")

after_condition_umap <- DimPlot(final_obj, reduction = "umap.harmony",
                                 group.by = "condition") +
  ggtitle("iPSC_PD: After Harmony Integration (Condition)") +
  labs(color = "Condition")

timepoint_umap <- DimPlot(final_obj, reduction = "umap.harmony",
                          group.by = "Timepoint", pt.size = 1.5) +
  ggtitle("iPSC_PD: Timepoint (DIV stage)") +
  labs(color = "Timepoint") +
  scale_color_viridis_d(option = "magma")

pca_donor_plot <- DimPlot(final_obj, reduction = "pca",
                          dims = c(1, 2),
                          group.by = DONOR_COL) +
  ggtitle("iPSC_PD: PCA Colored by Donor (pre-Harmony)")

# =============================================================================
# STEP 8: Generate 7-page PDF report
# =============================================================================
cat("\n=== Generating PDF report:", OUTPUT_PDF, "===\n")

pdf(OUTPUT_PDF, width = 14, height = 10)

# Page 1: PCA Elbow Analysis
grid.arrange(
  pre_integration$elbow_plots$elbow_plot,
  pre_integration$elbow_plots$quant_plot,
  nrow = 1, ncol = 2,
  top = textGrob("iPSC_PD — PCA Elbow Analysis (Pre-Integration)",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Page 2: Before / After Harmony (by donor and condition, 2x2)
grid.arrange(
  harmony_result$plots$before,
  harmony_result$plots$after,
  before_condition_umap,
  after_condition_umap,
  nrow = 2, ncol = 2,
  top = textGrob("iPSC_PD — Harmony Integration Comparison (Donor & Condition)",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Page 3: Clustree
grid.arrange(
  clustering_result$plots$clustree,
  nrow = 1, ncol = 1,
  top = textGrob("iPSC_PD — Clustering Resolution Tree",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Page 4: Final clusters + silhouette table
grid.arrange(
  clustering_result$plots$cluster_umap,
  tableGrob(head(clustering_result$sil_results, 10)),
  nrow = 1, ncol = 2,
  top = textGrob("iPSC_PD — Final Clustering Results",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Page 5: QC metrics (2x2 grid)
grid.arrange(
  qc_plots$cellcycle,
  qc_plots$nfeature,
  qc_plots$ncount,
  clustering_result$plots$cluster_umap,
  nrow = 2, ncol = 2,
  top = textGrob("iPSC_PD — QC Metrics on Harmony UMAP",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Page 6: Donor + condition UMAPs (batch diagnostics)
grid.arrange(
  donor_umap,
  condition_umap,
  nrow = 1, ncol = 2,
  top = textGrob("iPSC_PD — Donor and Condition Distribution",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Page 7: Timepoint (DIV stage) UMAP
grid.arrange(
  timepoint_umap,
  nrow = 1, ncol = 1,
  top = textGrob("iPSC_PD — Timepoint (DIV Stage) Distribution",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

dev.off()
cat("  PDF saved:", OUTPUT_PDF, "\n")

# =============================================================================
# STEP 9: Save integrated object + summary
# =============================================================================
cat("\n=== Saving integrated Seurat object ===\n")
saveRDS(final_obj, file = OUTPUT_RDS)
cat("  Saved:", OUTPUT_RDS, "\n")

cat(sprintf("
========================================
INTEGRATION COMPLETE — iPSC_PD
========================================
  Input           : %s
  Output RDS      : %s
  Output PDF      : %s

  Total cells     : %d
  Optimal PCs     : %d
  Calc. optimal resolution : %s
  Final resolution used    : %s (%s)
  Number of clusters       : %d

  Cells per condition:
", INPUT_RDS, OUTPUT_RDS, OUTPUT_PDF,
   ncol(final_obj),
   pre_integration$optimal_pcs,
   clustering_result$optimal_res,
   clustering_result$final_res,
   clustering_result$resolution_method,
   length(unique(final_obj$seurat_clusters))))

print(table(Condition = final_obj@meta.data$condition))
cat("\n  Cells per donor:\n")
print(table(Donor = final_obj@meta.data[[DONOR_COL]]))
