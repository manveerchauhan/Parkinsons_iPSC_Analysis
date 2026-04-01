#!/usr/bin/env Rscript
# =============================================================================
# qc_filter.R
# =============================================================================
# Per-donor adaptive QC filtering for PPMI iPSC scRNA-seq data.
#
# Adapted from: ch3_org_working_analysis.R (Neurodevelopmental_Models_analysis)
# Author : Manveer Chauhan
# Usage  : Rscript qc_filter.R  (run from Parkinsons_iPSC_Analysis/)
#
# Input  : export/seurat_filtered.rds
# Output : export/seurat_qc_filtered.rds
#          output_files/qc/iPSC_QC_filter_report.pdf
#
# Strategy: Split merged object by donor → adaptive threshold per donor →
#           merge filtered slices back. Thresholds are derived from density
#           peak detection (rightmost of top-2 peaks by height) ± 1.5×SD.
#           All thresholds can be overridden per-donor in the qc_params table.
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_RDS  <- "export/seurat_filtered.rds"
OUTPUT_RDS <- "export/seurat_qc_filtered.rds"
OUTPUT_PDF <- "output_files/qc/iPSC_QC_filter_report.pdf"
DONOR_COL  <- "PPMI Donor ID"

# Per-donor override table — NA = use adaptive threshold derived from density peaks.
# After reviewing qc_visualisation.R output, edit specific donor rows to override.
# Columns:
#   custom_min_genes  : override MIN_GENES (NA = adaptive: peak - 1.5*SD)
#   custom_max_genes  : override MAX_GENES (NA = adaptive: peak + 1.5*SD)
#   custom_max_counts : override MAX_COUNTS (NA = adaptive: peak + 1.5*SD)
#   min_counts        : hard floor for nCount_RNA (200 = conservative for pre-QC'd data)
#   mt_threshold      : mitochondrial % upper limit (fixed at 7.5 for all donors)
DONOR_OVERRIDES <- tibble::tibble(
  donor             = character(0),
  custom_min_genes  = double(0),
  custom_max_genes  = double(0),
  custom_max_counts = double(0),
  min_counts        = double(0),
  mt_threshold      = double(0)
)
# To add a custom override, uncomment and edit a row like this:
# DONOR_OVERRIDES <- tibble::add_row(DONOR_OVERRIDES,
#   donor = "3406", custom_min_genes = NA, custom_max_genes = 8000,
#   custom_max_counts = NA, min_counts = 200, mt_threshold = 7.5)

# =============================================================================
# LIBRARIES
# =============================================================================
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)

# =============================================================================
# HELPERS
# =============================================================================
section <- function(title) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  ", title, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

# =============================================================================
# QC FUNCTIONS (adapted from ch3_org_working_analysis.R)
# =============================================================================

#' Detect dominant peak in a numeric distribution.
#' Finds all local maxima, takes the top-2 by density height,
#' then returns the rightmost (highest x-value) of those two.
find_dominant_peak <- function(values, bw = "nrd0") {
  dens <- density(values, bw = bw, na.rm = TRUE)

  peak_idx <- which(diff(sign(diff(dens$y))) == -2) + 1

  if (length(peak_idx) == 0) {
    return(list(
      center       = mean(values, na.rm = TRUE),
      method       = "mean",
      n_peaks_total = 0,
      n_peaks_top2  = 0,
      all_peaks    = NULL,
      top2_peaks   = NULL,
      selected_peak = mean(values, na.rm = TRUE),
      density_obj  = dens
    ))
  }

  all_peak_values    <- dens$x[peak_idx]
  all_peak_densities <- dens$y[peak_idx]

  if (length(peak_idx) == 1) {
    peak_value       <- all_peak_values[1]
    top2_peak_values <- all_peak_values[1]
    n_top2           <- 1
  } else {
    top2_idx         <- order(all_peak_densities, decreasing = TRUE)[1:min(2, length(all_peak_densities))]
    top2_peak_values <- all_peak_values[top2_idx]
    n_top2           <- length(top2_peak_values)
    peak_value       <- max(top2_peak_values)
  }

  return(list(
    center        = peak_value,
    method        = "peak",
    n_peaks_total = length(peak_idx),
    n_peaks_top2  = n_top2,
    all_peaks     = all_peak_values,
    top2_peaks    = top2_peak_values,
    selected_peak = peak_value,
    density_obj   = dens
  ))
}

#' Plot distribution histogram with density overlay and threshold lines.
plot_threshold_diagnostics <- function(values, metric_name, peak_result,
                                       min_threshold, max_threshold, donor_id) {
  plot_df <- data.frame(value = values)

  p <- ggplot(plot_df, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = "lightblue", alpha = 0.6) +
    geom_density(color = "black", linewidth = 1) +
    theme_classic() +
    labs(
      title    = paste0(donor_id, ": ", metric_name, " Distribution"),
      subtitle = paste0("Method: ", peak_result$method,
                        " | Selected peak: ", round(peak_result$center, 1),
                        " | Total peaks: ", peak_result$n_peaks_total,
                        " | Top 2 peaks used"),
      x = metric_name,
      y = "Density"
    )

  if (!is.null(peak_result$all_peaks) && length(peak_result$all_peaks) > 0) {
    all_dens <- approx(peak_result$density_obj$x, peak_result$density_obj$y,
                       xout = peak_result$all_peaks)$y
    p <- p + geom_point(data = data.frame(x = peak_result$all_peaks, y = all_dens),
                        aes(x = x, y = y), color = "gray60", size = 2, shape = 19)
  }

  if (!is.null(peak_result$top2_peaks) && length(peak_result$top2_peaks) > 0) {
    top2_dens <- approx(peak_result$density_obj$x, peak_result$density_obj$y,
                        xout = peak_result$top2_peaks)$y
    p <- p + geom_point(data = data.frame(x = peak_result$top2_peaks, y = top2_dens),
                        aes(x = x, y = y), color = "orange", size = 4, shape = 19)
  }

  y_max <- max(peak_result$density_obj$y) * 0.95
  p <- p +
    geom_vline(xintercept = peak_result$center,
               color = "darkgreen", linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = min_threshold,
               color = "red", linetype = "dotted", linewidth = 0.8) +
    geom_vline(xintercept = max_threshold,
               color = "red", linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = min_threshold, y = y_max,
             label = paste0("MIN: ", round(min_threshold)),
             angle = 90, vjust = -0.5, size = 3, color = "red") +
    annotate("text", x = max_threshold, y = y_max,
             label = paste0("MAX: ", round(max_threshold)),
             angle = 90, vjust = 1.5, size = 3, color = "red")

  return(p)
}

#' Apply adaptive per-donor QC filtering with before/after diagnostics.
#' @param seurat_donor  Single-donor Seurat object (from SplitObject)
#' @param donor_params  One-row tibble from qc_params for this donor
filter_donor <- function(seurat_donor, donor_params) {

  donor_id      <- donor_params$donor
  initial_cells <- ncol(seurat_donor)

  cat("\n  --- Donor:", donor_id, "---\n")
  cat("    Initial cells:", initial_cells, "\n")

  # Ensure percent.mt is present (^MT- for PPMI HGNC symbols)
  if (!"percent.mt" %in% colnames(seurat_donor@meta.data)) {
    seurat_donor[["percent.mt"]] <- PercentageFeatureSet(seurat_donor, pattern = "^MT-")
  }

  # ---- nFeature_RNA thresholds ----
  if (is.na(donor_params$custom_min_genes) || is.na(donor_params$custom_max_genes)) {
    peak_result_genes <- find_dominant_peak(seurat_donor$nFeature_RNA)
    center_genes      <- peak_result_genes$center
    sd_genes          <- sd(seurat_donor$nFeature_RNA)
    cat("    nFeature_RNA — peak:", round(center_genes),
        "| SD:", round(sd_genes),
        "| n_peaks:", peak_result_genes$n_peaks_total, "\n")
  }

  MIN_GENES <- if (!is.na(donor_params$custom_min_genes)) {
    cat("    MIN_GENES: CUSTOM =", donor_params$custom_min_genes, "\n")
    peak_result_genes <- NULL
    donor_params$custom_min_genes
  } else {
    val <- round(center_genes - 1.25 * sd_genes)
    cat("    MIN_GENES: adaptive =", val, "\n")
    val
  }

  MAX_GENES <- if (!is.na(donor_params$custom_max_genes)) {
    cat("    MAX_GENES: CUSTOM =", donor_params$custom_max_genes, "\n")
    donor_params$custom_max_genes
  } else {
    val <- round(center_genes + 1.25 * sd_genes)
    cat("    MAX_GENES: adaptive =", val, "\n")
    val
  }

  # ---- nCount_RNA thresholds ----
  if (is.na(donor_params$custom_max_counts)) {
    peak_result_counts <- find_dominant_peak(seurat_donor$nCount_RNA)
    center_counts      <- peak_result_counts$center
    sd_counts          <- sd(seurat_donor$nCount_RNA)
    cat("    nCount_RNA  — peak:", round(center_counts),
        "| SD:", round(sd_counts),
        "| n_peaks:", peak_result_counts$n_peaks_total, "\n")
    MAX_COUNTS <- round(center_counts + 1.25 * sd_counts)
    cat("    MAX_COUNTS: adaptive =", MAX_COUNTS, "\n")
  } else {
    peak_result_counts <- NULL
    MAX_COUNTS <- donor_params$custom_max_counts
    cat("    MAX_COUNTS: CUSTOM =", MAX_COUNTS, "\n")
  }

  MIN_COUNTS   <- donor_params$min_counts
  MT_THRESHOLD <- donor_params$mt_threshold
  cat("    MIN_COUNTS:", MIN_COUNTS, "| MT_THRESHOLD:", MT_THRESHOLD, "%\n")

  # ---- Before-filtering plots ----
  vln_before <- plot_grid_vlns(seurat_donor, donor_id, tag = "Before")
  scatter_before <- FeatureScatter(seurat_donor,
                                   feature1 = "nCount_RNA",
                                   feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm", formula = y ~ x) +
    ggtitle(paste0(donor_id, ": Before Filtering")) +
    NoLegend()

  # ---- Apply filter ----
  seurat_donor_filt <- subset(seurat_donor,
    subset = nFeature_RNA > MIN_GENES  &
             nFeature_RNA < MAX_GENES  &
             nCount_RNA   > MIN_COUNTS &
             nCount_RNA   < MAX_COUNTS &
             percent.mt   < MT_THRESHOLD)

  final_cells   <- ncol(seurat_donor_filt)
  cells_removed <- initial_cells - final_cells
  pct_kept      <- round((final_cells / initial_cells) * 100, 1)
  cat("    After filtering:", final_cells,
      "cells (removed", cells_removed, "/", 100 - pct_kept, "% )\n")

  # ---- After-filtering plots ----
  vln_after <- plot_grid_vlns(seurat_donor_filt, donor_id, tag = "After")
  scatter_after <- FeatureScatter(seurat_donor_filt,
                                  feature1 = "nCount_RNA",
                                  feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm", formula = y ~ x) +
    ggtitle(paste0(donor_id, ": After Filtering")) +
    NoLegend()

  # ---- Diagnostic distribution plots ----
  diag_plots <- list()
  if (!is.null(peak_result_genes)) {
    diag_plots$genes_dist <- plot_threshold_diagnostics(
      seurat_donor$nFeature_RNA, "nFeature_RNA (Genes per Cell)",
      peak_result_genes, MIN_GENES, MAX_GENES, donor_id)
  }
  if (!is.null(peak_result_counts)) {
    diag_plots$counts_dist <- plot_threshold_diagnostics(
      seurat_donor$nCount_RNA, "nCount_RNA (UMI Counts per Cell)",
      peak_result_counts, MIN_COUNTS, MAX_COUNTS, donor_id)
  }

  return(list(
    seurat_obj = seurat_donor_filt,
    plots = list(
      vln_before   = vln_before,
      vln_after    = vln_after,
      scatter_before = scatter_before,
      scatter_after  = scatter_after,
      diagnostics    = diag_plots
    ),
    thresholds = list(
      MIN_GENES    = MIN_GENES,
      MAX_GENES    = MAX_GENES,
      MIN_COUNTS   = MIN_COUNTS,
      MAX_COUNTS   = MAX_COUNTS,
      MT_THRESHOLD = MT_THRESHOLD
    ),
    summary = data.frame(
      donor         = donor_id,
      cells_before  = initial_cells,
      cells_after   = final_cells,
      cells_removed = cells_removed,
      pct_kept      = pct_kept,
      MIN_GENES     = MIN_GENES,
      MAX_GENES     = MAX_GENES,
      MIN_COUNTS    = MIN_COUNTS,
      MAX_COUNTS    = MAX_COUNTS,
      MT_THRESHOLD  = MT_THRESHOLD
    )
  ))
}

# Helper: 3-panel violin plot (nFeature / nCount / percent.mt) for one donor
plot_grid_vlns <- function(obj, donor_id, tag) {
  v1 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0.1) +
    ggtitle(paste0(donor_id, ": Genes — ", tag)) + NoLegend()
  v2 <- VlnPlot(obj, features = "nCount_RNA",   pt.size = 0.1) +
    ggtitle(paste0(donor_id, ": Counts — ", tag)) + NoLegend()
  v3 <- VlnPlot(obj, features = "percent.mt",   pt.size = 0.1) +
    ggtitle(paste0(donor_id, ": MT% — ", tag)) + NoLegend()
  list(v1, v2, v3)
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

# Add percent.mt to the full object (used for violin plots; reused inside filter_donor)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
cat("  percent.mt added (genes matching ^MT-)\n")

# =============================================================================
# BUILD qc_params TABLE
# =============================================================================
section("BUILDING QC PARAMETER TABLE")

all_donors <- sort(unique(seurat_obj@meta.data[[DONOR_COL]]))
cat("  Donors found:", paste(all_donors, collapse = ", "), "\n")

# Default params for all donors (all adaptive)
qc_params <- tibble::tibble(
  donor             = all_donors,
  custom_min_genes  = NA_real_,
  custom_max_genes  = NA_real_,
  custom_max_counts = NA_real_,
  min_counts        = 200,
  mt_threshold      = 7.5
)

# Apply any user-specified overrides
if (nrow(DONOR_OVERRIDES) > 0) {
  for (i in seq_len(nrow(DONOR_OVERRIDES))) {
    ov  <- DONOR_OVERRIDES[i, ]
    idx <- which(qc_params$donor == ov$donor)
    if (length(idx) == 0) {
      warning("Override donor '", ov$donor, "' not found in data — skipping.")
      next
    }
    qc_params[idx, ] <- ov
    cat("  Applied custom override for donor:", ov$donor, "\n")
  }
}

cat("\n  QC parameter table:\n")
print(as.data.frame(qc_params))

# =============================================================================
# SPLIT → FILTER PER DONOR → MERGE
# =============================================================================
section("PER-DONOR ADAPTIVE FILTERING")

donor_list <- SplitObject(seurat_obj, split.by = DONOR_COL)
cat("  Split into", length(donor_list), "donor objects\n")

filter_results <- lapply(names(donor_list), function(d) {
  params <- qc_params %>% dplyr::filter(donor == d)
  filter_donor(donor_list[[d]], params)
})
names(filter_results) <- names(donor_list)

# Merge filtered objects back
cat("\n  Merging filtered donor objects...\n")
filtered_objs  <- lapply(filter_results, function(r) r$seurat_obj)
seurat_qc_filt <- merge(filtered_objs[[1]],
                         y          = filtered_objs[-1],
                         add.cell.ids = names(filtered_objs))
cat("  Merged object cells:", ncol(seurat_qc_filt), "\n")

# Seurat merge() applies make.names() to metadata columns (spaces → dots).
# Restore the original column name so all downstream scripts work unchanged.
safe_col <- make.names(DONOR_COL)   # "PPMI.Donor.ID"
if (safe_col %in% colnames(seurat_qc_filt@meta.data) &&
    !DONOR_COL %in% colnames(seurat_qc_filt@meta.data)) {
  idx <- which(colnames(seurat_qc_filt@meta.data) == safe_col)
  colnames(seurat_qc_filt@meta.data)[idx] <- DONOR_COL
  cat("  Restored metadata column:", DONOR_COL,
      "(was mangled to:", safe_col, "by merge)\n")
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================
section("FILTERING SUMMARY")

summary_df <- do.call(rbind, lapply(filter_results, function(r) r$summary))
summary_df <- summary_df %>%
  dplyr::left_join(
    seurat_obj@meta.data %>%
      dplyr::select(donor = all_of(DONOR_COL), condition) %>%
      dplyr::distinct(),
    by = "donor"
  ) %>%
  dplyr::select(donor, condition, everything()) %>%
  dplyr::arrange(condition, donor)

cat("\n")
print(as.data.frame(summary_df))

cat(sprintf("\n  Total: %d → %d cells (removed %d, %.1f%% kept)\n",
    ncol(seurat_obj),
    ncol(seurat_qc_filt),
    ncol(seurat_obj) - ncol(seurat_qc_filt),
    100 * ncol(seurat_qc_filt) / ncol(seurat_obj)))

# =============================================================================
# SAVE PDF REPORT
# =============================================================================
section("SAVING PDF REPORT")

if (!dir.exists(dirname(OUTPUT_PDF))) {
  dir.create(dirname(OUTPUT_PDF), recursive = TRUE)
}

pdf(OUTPUT_PDF, width = 14, height = 10)

for (d in names(filter_results)) {
  res <- filter_results[[d]]

  # Page A: before/after violin plots (3 before + 3 after + 2 scatters)
  grid.arrange(
    res$plots$vln_before[[1]], res$plots$vln_before[[2]], res$plots$vln_before[[3]],
    res$plots$vln_after[[1]],  res$plots$vln_after[[2]],  res$plots$vln_after[[3]],
    nrow = 2, ncol = 3,
    top = textGrob(paste0(d, " — QC Violin Plots (Before / After)"),
                   gp = gpar(fontsize = 13, fontface = "bold"))
  )

  # Page B: scatter before/after + diagnostic density plots
  diag_list <- res$plots$diagnostics
  plot_list  <- list(res$plots$scatter_before, res$plots$scatter_after)
  if (!is.null(diag_list$genes_dist))  plot_list <- c(plot_list, list(diag_list$genes_dist))
  if (!is.null(diag_list$counts_dist)) plot_list <- c(plot_list, list(diag_list$counts_dist))

  n_cols <- min(2, length(plot_list))
  n_rows <- ceiling(length(plot_list) / n_cols)

  grid.arrange(
    grobs = plot_list,
    nrow = n_rows, ncol = n_cols,
    top = textGrob(paste0(d, " — Scatter & Threshold Diagnostics"),
                   gp = gpar(fontsize = 13, fontface = "bold"))
  )
}

# Final summary page
grid.arrange(
  tableGrob(as.data.frame(summary_df),
            rows = NULL,
            theme = ttheme_default(base_size = 9)),
  nrow = 1,
  top = textGrob("iPSC QC Filtering Summary — All Donors",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

dev.off()
cat("  Saved:", OUTPUT_PDF, "\n")

# =============================================================================
# SAVE FILTERED SEURAT OBJECT
# =============================================================================
section("SAVING FILTERED SEURAT OBJECT")

saveRDS(seurat_qc_filt, file = OUTPUT_RDS)
cat("  Saved:", OUTPUT_RDS, "\n")

cat(sprintf("
======================================================================
QC FILTERING COMPLETE
======================================================================
  Input           : %s  (%d cells)
  Output RDS      : %s  (%d cells)
  Output PDF      : %s

  Cells removed   : %d  (%.1f%% of input)
  Donors filtered : %d

  Review %s to check:
    - Threshold lines align with distribution peaks
    - No donor loses >50%% of cells (check pct_kept column)
    - Donors with unusual distributions may need custom overrides in
      the DONOR_OVERRIDES table at the top of this script.

",
INPUT_RDS,  ncol(seurat_obj),
OUTPUT_RDS, ncol(seurat_qc_filt),
OUTPUT_PDF,
ncol(seurat_obj) - ncol(seurat_qc_filt),
100 * (1 - ncol(seurat_qc_filt) / ncol(seurat_obj)),
length(donor_list),
OUTPUT_PDF))
