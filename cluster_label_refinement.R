# ==============================================================================
# Cluster Label Refinement — iPSC PD Integrated Harmony Object
# ==============================================================================
# Author: Manveer Chauhan
# 
# Purpose: Replace scType automatic annotations with manually curated labels
#          based on systematic DotPlot analysis of marker gene panels.
#
# Rationale: scType's default brain reference database is poorly suited for 
#            iPSC-derived midbrain cultures (Kriks et al. 2011 protocol).
#            Labels such as "Cancer cells" and "Non myelinating Schwann cells"
#            are artefacts — driven by proliferation genes and peripheral nervous
#            system markers in the reference, respectively.
#
# Evidence: Each revised label below is supported by expression patterns across
#           9 marker panels (Pan-Neuronal, Immature Neuron, Mature Neuron,
#           Midbrain DA Lineage, Progenitor/Stem Cell, Forebrain Off-Target,
#           Hypothalamic Off-Target, Neurotransmitter Subtypes, and Glial).
#
# Input  : output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds
# Output : output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds (overwritten)
#          output_files/cell_annotation/annotation_progression_UMAP.pdf
#          output_files/cell_annotation/refined_labels_UMAP_full.pdf
# ==============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("/data/gpfs/projects/punim2251/Parkinsons_iPSC_Analysis")
set.seed(5728)

dir.create("./output_files/cell_annotation", recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("CLUSTER LABEL REFINEMENT\n")
cat("========================================\n\n")

# Load Harmony-integrated object
seurat_obj <- readRDS("output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds")
cat("Loaded object:", ncol(seurat_obj), "cells,", nlevels(seurat_obj$seurat_clusters), "clusters\n\n")

# ==============================================================================
# DEFINE REVISED LABELS
# ==============================================================================
# Mapping: scType label -> revised annotation
#
# Provenance: sctype_ScTypeDB_Default (original, untouched)
#               → cell_type_refined (copy, then remapped)
#
# Each revised label is grounded in marker gene DotPlot evidence interpreted
# in the context of the Kriks et al. (2011) floor-plate-based midbrain
# dopaminergic neuron differentiation protocol.
#
# Evidence summary (by scType label):
#
#  "Non myelinating Schwann cells" → "Immature midbrain neurons"  [cluster 0]
#    Pan-neuronal+++, DCX++, SLC17A6+, LMX1A weak — neuronal, not Schwann.
#    Progenitors largely exited (SOX2 low). Not yet synaptically mature.
#    No forebrain markers. Artefact from peripheral NS genes in scType DB.
#
#  "Cancer cells" → "Neural progenitors"  [cluster 1]
#    SOX2+++, VIM+, NES+ — clear progenitor signature.
#    Pan-neuronal near-absent; immature neuron markers near zero.
#    "Cancer" call driven by proliferation genes in scType DB.
#
#  "Unknown" → "Immature neurons (transitioning)"  [cluster 2]
#    Pan-neuronal++, DCX++; progenitor markers downregulating.
#    Profile similar to cluster 0 but slightly less mature.
#
#  "Mature neurons" → "Mature neurons (glutamatergic)"  [cluster 4]
#    RBFOX3 strongest of any cluster; SLC17A6++ confirms glutamatergic character.
#    KCNJ6 moderate; DA markers weak. Refined label adds neurotransmitter identity.
#
#  "Dopaminergic neurons" → "Midbrain dopaminergic neurons"  [cluster 7]
#    FOXA2 strongest in dataset, NR4A2+++, TH++, KCNJ6+ — midbrain DA.
#    SLC6A3 near-absent (DAT onset ~3-4 months; consistent with Kriks et al.).
#    Refined label adds regional qualifier (midbrain, not hypothalamic DA).
#
#  "Radial glial cells" → "Midbrain floor plate progenitors"  [clusters 5, 9, 10]
#    SOX2+++, VIM+++, NES+++, OTX2+++ (strongest anywhere), LMX1A++.
#    Defines the FOXA2+/LMX1A+/OTX2+ FP precursor of the Kriks protocol.
#    Cluster 10 is small/isolated; shares same progenitor profile.
#    TODO: run FindMarkers(ident.1=10, ident.2=c(5,9)) to confirm.
#
#  "GABAergic neurons" → "Off-target forebrain GABAergic neurons"  [cluster 6]
#    GAD1/2 strongly and exclusively expressed. EMX2+++, LHX2+++, SIX3++.
#    FOXA2 absent — not midbrain-derived. Predicted by Kriks et al. Fig. 2d
#    as LSB/S/F8 condition output (incomplete CHIR patterning).
#
#  "Astrocytes" → "Astrocytes" (confirmed)  [cluster 8]
#    GFAP and S100B exclusively expressed; only cluster with glial markers.
#
#  "Glutamatergic neurons" → "Mature glutamatergic neurons"  [clusters 3, 7*]
#    *Note: scType assigned "Glutamatergic neurons" to clusters 3 AND 7.
#    Cluster 7 is the DA cluster (see Dopaminergic neurons above); scType
#    misclassified it. Both clusters sharing this label receive the same
#    refined label here (Option A) because cluster 7 is disambiguated by
#    the "Dopaminergic neurons" remap above, which runs first.
#    Cluster 3: SNAP25++, SYP++, RBFOX3++, SLC17A6+++ — confirmed mature glut.

sctype_to_refined <- c(
  "Non myelinating Schwann cells" = "Immature midbrain neurons",
  "Unknown"                       = "Immature neurons (transitioning)",
  "Mature neurons"                = "Mature neurons (glutamatergic)",
  "Dopaminergic neurons"          = "Midbrain dopaminergic neurons",
  "GABAergic neurons"             = "Off-target forebrain GABAergic neurons",
  "Astrocytes"                    = "Astrocytes",
  "Cancer cells"                  = "Neural progenitors",
  "Radial glial cells"            = "Midbrain floor plate progenitors"
)
# "Glutamatergic neurons" is handled separately after the loop (see below)

# ==============================================================================
# APPLY REVISED LABELS
# ==============================================================================

cat("Applying revised cluster labels...\n\n")

# Store the original scType labels for reference (do not overwrite)
if (!"sctype_original" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj@meta.data$sctype_original <- seurat_obj@meta.data$sctype_ScTypeDB_Default
  cat("  Original scType labels preserved in 'sctype_original' column\n")
}

# Duplicate scType column as starting point — refined labels derived from it
seurat_obj@meta.data$cell_type_refined <- as.character(seurat_obj@meta.data$sctype_ScTypeDB_Default)

# Remap unambiguous scType labels (1 scType label → 1 revised label)
for (old_label in names(sctype_to_refined)) {
  mask <- seurat_obj@meta.data$cell_type_refined == old_label
  seurat_obj@meta.data$cell_type_refined[mask] <- sctype_to_refined[[old_label]]
  cat("  Remapped:", old_label, "→", sctype_to_refined[[old_label]],
      paste0("(", sum(mask), " cells)\n"))
}

# Handle "Glutamatergic neurons" — shared by clusters 3 & 7; both → same refined label
glut_mask <- seurat_obj@meta.data$cell_type_refined == "Glutamatergic neurons"
seurat_obj@meta.data$cell_type_refined[glut_mask] <- "Mature glutamatergic neurons"
cat("  Remapped: Glutamatergic neurons → Mature glutamatergic neurons",
    paste0("(", sum(glut_mask), " cells)\n"))

# Set factor order: progenitors → immature → mature → off-target → glia
label_order <- c(
  "Neural progenitors",
  "Midbrain floor plate progenitors",
  "Immature neurons (transitioning)",
  "Immature midbrain neurons",
  "Mature neurons (glutamatergic)",
  "Mature glutamatergic neurons",
  "Midbrain dopaminergic neurons",
  "Off-target forebrain GABAergic neurons",
  "Astrocytes"
)
seurat_obj@meta.data$cell_type_refined <- factor(seurat_obj@meta.data$cell_type_refined,
                                                  levels = label_order)

# Set active Idents to refined labels for downstream analysis
Idents(seurat_obj) <- "cell_type_refined"

cat("\n  Cells per refined label:\n")
print(table(seurat_obj@meta.data$cell_type_refined))

# Verification table — 3-level mapping: cluster → scType → refined
cat("\n  Verification (cluster → scType → refined):\n")
verification <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, sctype_ScTypeDB_Default, cell_type_refined) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(as.numeric(as.character(seurat_clusters)))
print(as.data.frame(verification), row.names = FALSE)

# ==============================================================================
# GENERATE UMAP VISUALISATIONS
# ==============================================================================

cat("\nGenerating UMAP plots...\n")

# Colour palette — designed for interpretability:
#   Blues/teals  = progenitors
#   Greens       = immature neurons
#   Warm colours = mature neurons
#   Red          = DA neurons (the target population)
#   Purple       = off-target
#   Grey         = glia
# refined_colours <- c(
#   "Neural progenitors"             = "gray",  # light blue
#   "Midbrain floor plate progenitors"         = "#1F78B4",  # dark blue
#   "Immature neurons (transitioning)"         = "#B2DF8A",  # light green
#   "Immature midbrain neurons"                = "#33A02C",  # dark green
#   "Mature neurons (glutamatergic)"           = "#FDBF6F",  # light orange
#   "Mature glutamatergic neurons"             = "#FF7F00",  # orange
#   "Midbrain dopaminergic neurons"            = "#E31A1C",  # red
#   "Off-target forebrain GABAergic neurons"   = "#CAB2D6",  # lavender
#   "Astrocytes"                               = "#A6CEE3"   # grey
# )

# Panel 1: Unsupervised clusters (numbers only)
p1 <- DimPlot(seurat_obj,
              reduction = "umap.harmony",
              group.by = "seurat_clusters",
              label = TRUE,
              label.size = 5,
              repel = TRUE) +
  labs(title = "Unsupervised Clusters",
       subtitle = "Harmony integration, resolution 0.2") +
  NoLegend() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9))

# Panel 2: Original scType labels
p2 <- DimPlot(seurat_obj,
              reduction = "umap.harmony",
              group.by = "sctype_ScTypeDB_Default",
              label = TRUE,
              label.size = 3,
              repel = TRUE) +
  labs(title = "scType Automatic Annotation",
       subtitle = "Reference: ScTypeDB_Default (Brain)") +
  NoLegend() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9))

# Panel 3: Refined labels
p3 <- DimPlot(seurat_obj,
              reduction = "umap.harmony",
              group.by = "cell_type_refined",
              label = TRUE,
              label.size = 3,
              repel = TRUE) +
              # cols = refined_colours) +
  labs(title = "Refined Annotations",
       subtitle = "Manual curation (Kriks et al. 2011 protocol context)") +
  NoLegend() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9))

# 3-panel annotation progression PDF
pdf("output_files/cell_annotation/annotation_progression_UMAP.pdf",
    width = 21, height = 7)
print(p1 + p2 + p3 +
        plot_annotation(
          title = "iPSC PD — Cell Type Annotation Progression",
          subtitle = "Cluster IDs → Automated (scType) → Manually Curated",
          theme = theme(
            plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40")
          )
        ))
dev.off()
cat("  Saved: output_files/cell_annotation/annotation_progression_UMAP.pdf\n")

# Full-size refined UMAP with legend
p_full <- DimPlot(seurat_obj,
                  reduction = "umap.harmony",
                  group.by = "cell_type_refined",
                  label = TRUE,
                  label.size = 3.5,
                  repel = TRUE) +
                  # cols = refined_colours) +
  labs(title = "iPSC PD: Refined Cell Type Annotations",
       subtitle = "Based on marker gene DotPlot analysis (Kriks et al. 2011 protocol context)",
       color = "Cell Type") +
  theme(plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.text   = element_text(size = 9))

# UMAP split by condition
p_split <- DimPlot(seurat_obj,
                   reduction = "umap.harmony",
                   group.by = "cell_type_refined",
                   split.by = "condition",
                   label = TRUE,
                   label.size = 3,
                   repel = TRUE) +
  labs(title = "iPSC PD: Refined Annotations by Condition",
       subtitle = "Split by condition metadata",
       color = "Cell Type") +
  theme(plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.text   = element_text(size = 9))

pdf("output_files/cell_annotation/refined_labels_UMAP_full.pdf",
    width = 18, height = 8)
print(p_full)
print(p_split)
dev.off()
cat("  Saved: output_files/cell_annotation/refined_labels_UMAP_full.pdf\n")

# ==============================================================================
# SAVE UPDATED OBJECT
# ==============================================================================

cat("\nSaving updated Seurat object...\n")
saveRDS(seurat_obj, "output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds")
cat("  Saved with new metadata column: 'cell_type_refined'\n")
cat("  Original scType labels preserved in: 'sctype_ScTypeDB_Default'\n")
cat("  Active Idents set to: 'cell_type_refined'\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("ANNOTATION REFINEMENT COMPLETE\n")
cat("========================================\n\n")

cat("Key changes from scType defaults:\n\n")
cat("  'Cancer cells'                  → Neural progenitors\n")
cat("  'Non myelinating Schwann cells' → Immature midbrain neurons\n")
cat("  'Unknown'                       → Immature neurons (transitioning)\n")
cat("  'Radial glial cells'            → Midbrain floor plate progenitors\n")
cat("  'GABAergic neurons'             → Off-target forebrain GABAergic neurons\n")
cat("  'Mature neurons'                → Mature neurons (glutamatergic)\n")
cat("  'Glutamatergic neurons'         → Mature glutamatergic neurons\n")
cat("  'Dopaminergic neurons'          → Midbrain dopaminergic neurons\n")
cat("  'Astrocytes'                    → Astrocytes (confirmed)\n\n")

cat("Output files:\n")
cat("  output_files/cell_annotation/annotation_progression_UMAP.pdf\n")
cat("  output_files/cell_annotation/refined_labels_UMAP_full.pdf (p1: all cells, p2: split by condition)\n")
cat("  output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds\n\n")

cat("Next steps:\n")
cat("  1. Verify cluster 10 assignment with FindMarkers(ident.1=10, ident.2=c(5,9))\n")
cat("  2. Run FindAllMarkers with refined labels for DE analysis\n")
cat("  3. Use 'cell_type_refined' as group.by for downstream comparisons\n")
cat("     between Idiopathic PD vs Control conditions\n")