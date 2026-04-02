## ============================================================
## iPSC_dotplot_markers.R
## Author: Manveer Chauhan
##
## Purpose: Generate a multi-page PDF of DotPlots for visual
##          verification of cell type annotations on the
##          Harmony-integrated iPSC Parkinson's Seurat object.
##
##          One DotPlot per biologically grouped marker gene set.
##          Clusters on x-axis, gene symbols on y-axis (coord_flip).
##
## Input:   output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds
## Output:  output_files/cell_annotation/marker_gene_dotplots.pdf
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

## ---- Paths --------------------------------------------------
INPUT_RDS  <- "output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds"
OUTPUT_PDF <- "output_files/cell_annotation/marker_gene_dotplots.pdf"

cat("\n####################################################\n")
cat("MARKER GENE DOTPLOT REPORT — iPSC PD\n")
cat("####################################################\n\n")

## ---- Load Seurat object ------------------------------------
cat("Loading annotated object:", INPUT_RDS, "\n")
if (!file.exists(INPUT_RDS)) {
  stop("Could not find: ", INPUT_RDS,
       "\n  Run iPSC_cell_annotation_sctype.R first.")
}
seurat_obj <- readRDS(INPUT_RDS)
cat("  Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "features\n")
cat("  Clusters:", length(unique(seurat_obj$seurat_clusters)), "\n\n")

## ---- Join layers (Seurat v5 compatibility) -----------------
cat("Joining layers for Seurat v5 compatibility...\n")
seurat_obj <- JoinLayers(seurat_obj)
cat("  Layers joined\n\n")

## ---- Set grouping column -----------------------------------
# Group by scType cell type labels (ScTypeDB_Default), falling back to
# seurat_clusters if the column is absent.
if ("sctype_ScTypeDB_Default" %in% colnames(seurat_obj@meta.data)) {
  cat("sctype_ScTypeDB_Default column found — grouping by cell type label\n")
  group_col <- "sctype_ScTypeDB_Default"
} else {
  cat("  WARNING: sctype_ScTypeDB_Default not found — using seurat_clusters\n")
  group_col <- "seurat_clusters"
}
group_col <- "cell_type_refined"
cat("  Grouping by:", group_col, "\n\n")

## ---- Gene symbol lookup ------------------------------------
# Features are stored as ENSG00000123456.10_SYMBOL or ENSG-SYMBOL.
# Strip everything before the last hyphen or underscore to get bare symbols.
cat("Building gene symbol lookup from feature names...\n")
all_features    <- rownames(seurat_obj)
feature_symbols <- sub(".*[-_]", "", all_features)
symbol_to_feature <- setNames(all_features, feature_symbols)
cat("  ", length(all_features), "features indexed\n\n")


## ---- Marker gene sets --------------------------------------
marker_sets <- list(

  "Pan-Neuronal Markers" = c(
    "TUBB3",   # beta-III-tubulin (TUJ1) — first neuronal marker at day 25 (Kriks et al. 2011)
    "MAP2",    # dendritic marker, increases with maturation
    "NCAM1",   # neural cell adhesion molecule — used for graft ID in Kriks et al.
    "STMN2"    # stathmin-2, robust pan-neuronal scRNA-seq marker (La Manno et al. 2016)
  ),

  "Immature Neuron Markers" = c(
    "DCX",     # doublecortin — upregulated day 13-25 in Kriks et al. Fig. 2c
    "ASCL1"    # MASH1 — upregulated day 13-25 in Kriks et al. Fig. 2c
  ),

  "Mature Neuron Markers" = c(
    "SYP",     # Synaptophysin — presynaptic vesicle marker
    "SNAP25",  # SNARE complex, synaptic vesicle fusion
    "MAPT",    # Tau — mature axonal marker, replaces DCX during maturation
    "RBFOX3"   # NeuN — mature neuronal nuclear marker
  ),

  "Midbrain DA Neuron Lineage (Kriks Protocol)" = c(
    "FOXA2",   # floor plate marker — THE key marker of the Kriks protocol (Fig. 1)
    "LMX1A",   # coexpressed with FOXA2 in midbrain FP (Fig. 1a,b)
    "OTX2",    # anterior/midbrain identity confirmation (Fig. 1a,c)
    "DDC",     # DOPA decarboxylase — present from day 11 (Fig. 1e)
    "TH",      # tyrosine hydroxylase — DA neuron marker from day 25 (Fig. 2a)
    "NR4A2",   # NURR1 — postmitotic DA marker (Fig. 2a,b)
    "PITX3",   # classic midbrain DA marker, expressed from day 25 (Fig. 2e)
    "SLC6A3",  # DAT — late maturation marker, onset 3-4 months (Fig. 4n)
    "KCNJ6",   # GIRK2/Kir3.2 — marks A9/SNpc-type DA neurons (Fig. 4o)
    "CALB1"    # Calbindin — marks A10/VTA-type DA neurons (Fig. 4p)
  ),

  "Progenitor and Stem Cell Markers" = c(
    "SOX2",    # pluripotency/NSC — should decrease along differentiation
    "VIM",     # vimentin — radial glia/neuroepithelial
    "NES"      # nestin — neural progenitor (used in Kriks Fig. 1a)
  ),

  "Off-Target: Forebrain Markers" = c(
    "FOXG1",   # CRITICAL — present in rosette-derived grafts, absent from FP-derived (Kriks et al.)
    "PAX6",    # dorsal forebrain progenitor (Kriks Fig. 1e,h)
    "EMX2",    # dorsal forebrain (Kriks Fig. 1e,h)
    "LHX2"     # dorsal forebrain (Kriks Fig. 1e,h)
  ),

  "Off-Target: Hypothalamic Markers" = c(
    "SIX3"     # hypothalamic (Kriks Fig. 1f)
  ),

  "Off-Target: Neurotransmitter Subtypes" = c(
    "GAD1",    # GABAergic — enriched in hypothalamic condition (Kriks Fig. 2d)
    "GAD2",    # GABAergic — complementary marker
    "SLC17A6", # VGLUT2 — expected in DA neurons as co-release marker (Hnasko et al. 2010)
    "SLC17A7", # VGLUT1 
    "SLC17A8", # VGLUT3
    "TPH2"     # serotonergic — Kriks found rare 5-HT+ neurons in FP cultures
  ),

  "Glial Markers" = c(
    "GFAP",    # astrocyte — mostly host-derived in Kriks grafts (Supp. Figs. 9-10)
    "S100B",    # astrocyte maturation marker
    "SLC1A3",
    "SLC1A2",
    "ALDH1L1",
    "SOX9"
  )
)

## ---- Generate PDF ------------------------------------------
cat("\nGenerating DotPlot PDF:", OUTPUT_PDF, "\n")

# Ensure output directory exists
dir.create(dirname(OUTPUT_PDF), showWarnings = FALSE, recursive = TRUE)

# DA lineage panel is larger (11 genes) so use a taller page for it
DA_SET_NAME <- "Midbrain DA Neuron Lineage (Kriks Protocol)"

pdf(OUTPUT_PDF, width = 14, height = 8)

for (set_name in names(marker_sets)) {

  cat("\n---", set_name, "---\n")

  symbols  <- marker_sets[[set_name]]
  resolved <- symbol_to_feature[symbols]
  found    <- resolved[!is.na(resolved)]
  missing  <- symbols[is.na(resolved)]

  cat("  Genes found:", length(found), "/", length(symbols), "\n")
  if (length(missing) > 0) {
    cat("  WARNING - genes not found in object:", paste(missing, collapse = ", "), "\n")
  }

  if (length(found) < 1) {
    cat("  SKIPPING — no genes found\n")
    next
  }

  # Use a taller page for the DA lineage panel
  if (set_name == DA_SET_NAME) {
    # Re-open with taller dimensions for this page only by using a blank page trick:
    # Seurat DotPlot renders to active device; we instead close/reopen the device
    # is not feasible mid-PDF. Instead adjust ggplot theme to compress vertically.
    # The 11-gene panel fits within height=8 using coord_flip (genes on y-axis).
    # If needed, reduce axis text size slightly for this panel.
    axis_text_size <- 10
  } else {
    axis_text_size <- 11
  }

  p <- DotPlot(
    seurat_obj,
    features = unname(found),
    group.by = group_col,
    dot.scale = 8,
    cols = c("lightgrey", "darkred")
  ) +
    # After coord_flip(), features are on the y-axis → strip ENSG prefix there
    scale_y_discrete(labels = function(x) sub(".*[-_]", "", x)) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y  = element_text(size = axis_text_size, face = "bold"),
      axis.title   = element_blank(),
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40"),
      legend.position = "right",
      plot.margin  = margin(t = 10, r = 10, b = 20, l = 10)
    ) +
    labs(
      title    = set_name,
      subtitle = "Dot size = % cells expressing  |  Colour = avg scaled expression"
    ) +
    coord_flip()

  print(p)
}

dev.off()

cat("\n####################################################\n")
cat("DONE\n")
cat("####################################################\n")
cat("PDF saved to:", OUTPUT_PDF, "\n\n")
