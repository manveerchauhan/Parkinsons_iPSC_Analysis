# ==============================================================================
# Pseudobulk Differential Expression — iPSC PD vs Control
# ==============================================================================
# Author: Manveer Chauhan
#
# Purpose: Identify differentially expressed genes between Idiopathic PD and
#          Control conditions using a pseudobulk DESeq2 approach. Cells are
#          aggregated per donor × condition × cell type before testing, avoiding
#          the inflated degrees-of-freedom problem of single-cell DE methods.
#
# Cell types tested:
#   - Astrocytes
#   - Midbrain dopaminergic neurons
#   - Mature glutamatergic neurons
#
# Input  : output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds
# Output : output_files/DEG/padj_histograms.pdf
# ==============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(DESeq2)

setwd("/data/gpfs/projects/punim2251/Parkinsons_iPSC_Analysis")

dir.create("output_files/DEG", recursive = TRUE, showWarnings = FALSE)

# Load Harmony-integrated object
seurat_obj <- readRDS("output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds")
cat("Loaded object:", ncol(seurat_obj), "cells,", nlevels(seurat_obj$seurat_clusters), "clusters\n\n")


PD.iPSC.pseudobulk <- AggregateExpression(seurat_obj, assays = "RNA",
                                       group.by = c("condition", "PPMI Donor ID", "cell_type_refined"),
                                       return.seurat = TRUE)

# ==============================================================================
# PSEUDOBULK DESeq2 — helper
# ==============================================================================
# For each cell type: subset pseudobulk samples → build DESeqDataSet → run DESeq2.
# This avoids the inflated-n problem of passing pseudobulk data through
# FindMarkers, which does not account for the donor structure.
#
# Filtering rule: keep genes with raw count >= 5 in at least 2 pseudobulk samples.
# This is intentionally lenient — independent filtering inside DESeq2 will
# further remove genes with insufficient mean expression.

run_pseudobulk_DESeq2 <- function(pb_obj, celltype, ref_condition = "Control") {

  cells_keep <- pb_obj$cell_type_refined == celltype
  sub        <- pb_obj[, cells_keep]

  cat(sprintf("\n--- %s ---\n", celltype))
  cat("Pseudobulk samples per condition:\n")
  print(table(sub$condition))

  cond_counts <- table(sub$condition)
  if (any(cond_counts < 2)) {
    warning(sprintf("Skipping '%s': fewer than 2 replicates in at least one condition.", celltype))
    return(NULL)
  }

  counts_mat <- GetAssayData(sub, assay = "RNA", layer = "counts")

  col_meta <- data.frame(
    condition = factor(sub$condition,
                       levels = c(ref_condition,
                                  setdiff(unique(sub$condition), ref_condition))),
    donor     = sub$`PPMI Donor ID`,
    row.names = colnames(sub)
  )

  dds <- DESeqDataSetFromMatrix(countData = round(counts_mat),
                                colData   = col_meta,
                                design    = ~ condition)

  keep <- rowSums(counts(dds) >= 5) >= 2
  cat(sprintf("Genes after pre-filtering: %d / %d\n", sum(keep), nrow(dds)))
  dds <- dds[keep, ]

  dds <- DESeq(dds)

  res <- results(dds,
                 contrast = c("condition", "Idiopathic-PD", ref_condition))
  cat(sprintf("Genes with padj < 0.05: %d\n", sum(res$padj < 0.05, na.rm = TRUE)))

  results_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)

  list(results = results_df, dds = dds)
}

# ==============================================================================
# RUN DEG — Idiopathic PD vs Control
# ==============================================================================

out.Astrocytes    <- run_pseudobulk_DESeq2(PD.iPSC.pseudobulk, "Astrocytes")
out.DA_neurons    <- run_pseudobulk_DESeq2(PD.iPSC.pseudobulk, "Midbrain dopaminergic neurons")
out.Gluta_neurons <- run_pseudobulk_DESeq2(PD.iPSC.pseudobulk, "Mature glutamatergic neurons")

Idiopathic.vs.CTRL.Astrocytes.pseudo    <- out.Astrocytes$results
Idiopathic.vs.CTRL.DA_neurons.pseudo    <- out.DA_neurons$results
Idiopathic.vs.CTRL.Gluta_neurons.pseudo <- out.Gluta_neurons$results

# ==============================================================================
# ADJUSTED P-VALUE HISTOGRAMS
# ==============================================================================

deg_results <- list(
  "Astrocytes"                      = Idiopathic.vs.CTRL.Astrocytes.pseudo,
  "Midbrain dopaminergic neurons"   = Idiopathic.vs.CTRL.DA_neurons.pseudo,
  "Mature glutamatergic neurons"    = Idiopathic.vs.CTRL.Gluta_neurons.pseudo
)

padj_plots <- lapply(names(deg_results), function(ct) {
  df <- deg_results[[ct]]
  if (is.null(df)) return(NULL)
  ggplot(df, aes(x = padj)) +
    geom_histogram(bins = 40, fill = "#4393C3", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    labs(title    = ct,
         subtitle = sprintf("%d genes tested, %d with padj < 0.05",
                            nrow(df), sum(df$padj < 0.05, na.rm = TRUE)),
         x        = "Adjusted p-value (BH)",
         y        = "Number of genes") +
    theme_classic() +
    theme(plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
})
padj_plots <- Filter(Negate(is.null), padj_plots)

pdf("output_files/DEG/padj_histograms.pdf", width = 14, height = 5)
print(wrap_plots(padj_plots, nrow = 1) +
        plot_annotation(title    = "iPSC PD — Pseudobulk DEG: Adjusted P-value Distributions",
                        subtitle = "Idiopathic PD vs Control | dashed line = padj 0.05",
                        theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                                         plot.subtitle = element_text(size = 10, colour = "grey40"))))
dev.off()
cat("Saved: output_files/DEG/padj_histograms.pdf\n")

# ==============================================================================
# DIAGNOSTICS — raw p-value distributions + dispersion plots
# ==============================================================================

dds_list <- list(
  "Astrocytes"                    = out.Astrocytes$dds,
  "Midbrain dopaminergic neurons" = out.DA_neurons$dds,
  "Mature glutamatergic neurons"  = out.Gluta_neurons$dds
)

pdf("output_files/DEG/DESeq2_diagnostics.pdf", width = 14, height = 5)

# Page 1: raw p-value histograms
pval_plots <- lapply(names(deg_results), function(ct) {
  df <- deg_results[[ct]]
  if (is.null(df)) return(NULL)
  ggplot(df, aes(x = pvalue)) +
    geom_histogram(bins = 40, fill = "#74C476", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    labs(title    = ct,
         subtitle = sprintf("%d genes, %d with raw p < 0.05",
                            sum(!is.na(df$pvalue)),
                            sum(df$pvalue < 0.05, na.rm = TRUE)),
         x = "Raw p-value",
         y = "Number of genes") +
    theme_classic() +
    theme(plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
})
pval_plots <- Filter(Negate(is.null), pval_plots)
print(wrap_plots(pval_plots, nrow = 1) +
        plot_annotation(title    = "iPSC PD — Raw P-value Distributions",
                        subtitle = "Enrichment near 0 indicates real signal; flat = null",
                        theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                                         plot.subtitle = element_text(size = 10, colour = "grey40"))))

# Pages 2–4: dispersion plots (one per cell type)
for (ct in names(dds_list)) {
  dds <- dds_list[[ct]]
  if (is.null(dds)) next
  plotDispEsts(dds, main = sprintf("%s — Dispersion Estimates", ct))
}

dev.off()
cat("Saved: output_files/DEG/DESeq2_diagnostics.pdf\n")
