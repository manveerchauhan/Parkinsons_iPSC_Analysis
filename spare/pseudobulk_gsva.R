#!/usr/bin/env Rscript
# =============================================================================
# pseudobulk_gsva.R
# =============================================================================
# Part 2: Pseudobulk ssGSEA analysis using the GSVA R package.
#
# Author : Manveer Chauhan
# Usage  : Rscript pseudobulk_gsva.R
#
# WORKFLOW:
#   A. Install / load libraries
#   B. Load scRNA-seq data into Seurat
#   C. Pseudobulk aggregation (AggregateExpression)
#   D. Normalize with edgeR (TMM → log-CPM)
#   E. Define gene sets + coverage check
#   F. Run ssGSEA (GSVA ≥ 2.0 with fallback to older syntax)
#   G. Visualizations (6 PDFs)
#   H. Statistical testing with limma
#   I. Final summary + save outputs
#
# PREREQUISITE: Run explore_and_export.py first, then fill in the TODO block.
# =============================================================================


# ============================================================
# TODO: Fill these in after running Part 1 (explore_and_export.py)
# ============================================================
H5AD_PATH     <- "export/adata_rawcounts.h5ad"  # path to exported h5ad
SAMPLE_COL    <- "TODO"    # column identifying biological replicates / donors
CELLTYPE_COL  <- "TODO"    # column with cell type annotations
CONDITION_COL <- "TODO"    # condition column (e.g. "disease" vs "control"), or NULL if none
                            # Example: CONDITION_COL <- NULL
COUNTS_LAYER  <- NULL      # set to layer name string if raw counts are in a layer,
                            # NULL if raw counts are already in X
                            # Example: COUNTS_LAYER <- "counts"

# Pseudobulk filtering
MIN_CELLS     <- 10        # minimum cells per pseudobulk group
MIN_COUNTS    <- 1000      # minimum total counts per pseudobulk group

# GSVA parameters
GSVA_METHOD   <- "ssgsea"  # "ssgsea" or "gsva"
NCORES        <- 4


# =============================================================================
# SECTION A: INSTALL AND LOAD LIBRARIES
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION A: Installing and loading libraries\n")
cat(strrep("=", 70), "\n\n")

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c("ggplot2", "pheatmap", "RColorBrewer", "tidyr", "dplyr",
               "Matrix", "reshape2", "cowplot")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing CRAN package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# ── Bioconductor packages ─────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c("GSVA", "edgeR", "BiocParallel", "limma")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing Bioconductor package:", pkg, "\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# ── Seurat ecosystem ──────────────────────────────────────────────────────────
seurat_pkgs <- c("Seurat", "anndataR", "SeuratDisk")
for (pkg in seurat_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    tryCatch({
      if (pkg == "SeuratDisk") {
        if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
        remotes::install_github("mojaveazure/seurat-disk", quiet = TRUE)
      } else if (pkg == "anndataR") {
        install.packages("anndataR", repos = "https://cloud.r-project.org")
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }, error = function(e) {
      cat("  WARNING: Could not install", pkg, "—", conditionMessage(e), "\n")
    })
  }
}

# ── Load all libraries ────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
  library(edgeR)
  library(limma)
  library(BiocParallel)
  library(GSVA)
})

cat("  All required libraries loaded.\n")
dir.create("gsva_results", showWarnings = FALSE)
cat("  Output directory 'gsva_results/' ready.\n")


# =============================================================================
# SECTION B: LOAD DATA INTO SEURAT
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION B: Loading data into Seurat\n")
cat(strrep("=", 70), "\n\n")

seurat_obj <- NULL

# ── Attempt 1: anndataR ───────────────────────────────────────────────────────
if (is.null(seurat_obj) && requireNamespace("anndataR", quietly = TRUE)) {
  cat("  Attempt 1: anndataR::read_h5ad() -> to_Seurat()\n")
  tryCatch({
    library(anndataR)
    adata <- read_h5ad(H5AD_PATH)
    seurat_obj <- adata$to_Seurat()
    cat("  SUCCESS with anndataR.\n")
  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
  })
}

# ── Attempt 2: SeuratDisk ─────────────────────────────────────────────────────
if (is.null(seurat_obj) && requireNamespace("SeuratDisk", quietly = TRUE)) {
  cat("  Attempt 2: SeuratDisk Convert + LoadH5Seurat\n")
  tryCatch({
    library(SeuratDisk)
    h5seurat_path <- sub("\\.h5ad$", ".h5seurat", H5AD_PATH)
    SeuratDisk::Convert(H5AD_PATH, dest = "h5seurat", overwrite = TRUE)
    seurat_obj <- SeuratDisk::LoadH5Seurat(h5seurat_path)
    cat("  SUCCESS with SeuratDisk.\n")
  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
  })
}

# ── Attempt 3: Manual MTX construction ───────────────────────────────────────
if (is.null(seurat_obj)) {
  cat("  Attempt 3: Manual construction from export/matrix.mtx\n")
  tryCatch({
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat is not installed.")
    }
    library(Seurat)

    cat("  Reading matrix.mtx ...\n")
    counts_matrix <- Matrix::readMM("export/matrix.mtx")  # genes × cells

    barcodes <- read.table("export/barcodes.tsv", header = FALSE,
                           stringsAsFactors = FALSE)[[1]]
    genes    <- read.table("export/genes.tsv", header = FALSE,
                           stringsAsFactors = FALSE)[[1]]
    metadata <- read.csv("export/metadata.csv", row.names = 1)

    rownames(counts_matrix) <- genes
    colnames(counts_matrix) <- barcodes

    seurat_obj <- CreateSeuratObject(counts  = counts_matrix,
                                     meta.data = metadata)
    cat("  SUCCESS with manual MTX construction.\n")
  }, error = function(e) {
    stop("All three Seurat loading methods failed: ", conditionMessage(e))
  })
}

# ── Handle COUNTS_LAYER ───────────────────────────────────────────────────────
if (!is.null(COUNTS_LAYER)) {
  cat("  Moving layer '", COUNTS_LAYER, "' to counts slot.\n", sep = "")
  tryCatch({
    library(Seurat)
    # Seurat v5 syntax
    LayerData(seurat_obj, layer = "counts") <-
      LayerData(seurat_obj, layer = COUNTS_LAYER)
  }, error = function(e) {
    tryCatch({
      # Seurat v4 syntax
      seurat_obj@assays$RNA@counts <-
        seurat_obj@assays$RNA@layers[[COUNTS_LAYER]]
    }, error = function(e2) {
      cat("  WARNING: Could not move layer to counts slot:", conditionMessage(e2), "\n")
    })
  })
}

# ── Diagnostics ───────────────────────────────────────────────────────────────
cat("\n  --- Seurat object diagnostics ---\n")
cat("  Dimensions:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n")
cat("  Assays:", paste(Assays(seurat_obj), collapse = ", "), "\n")

cat("\n  meta.data head:\n")
print(head(seurat_obj@meta.data, 3))

# Verify counts slot
tryCatch({
  cnt_check <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  cat("\n  Counts slot stats:\n")
  cnt_vals <- as.vector(cnt_check@x[seq_len(min(10000, length(cnt_check@x)))])
  cat("    min:", min(cnt_vals), "\n")
  cat("    max:", max(cnt_vals), "\n")
  cat("    is_integer:", all(cnt_vals == floor(cnt_vals)), "\n")
}, error = function(e) {
  cat("  WARNING: Could not inspect counts slot:", conditionMessage(e), "\n")
})

# Print unique values of key columns
for (col in c(SAMPLE_COL, CELLTYPE_COL, CONDITION_COL)) {
  if (!is.null(col) && col != "TODO" && col %in% colnames(seurat_obj@meta.data)) {
    vals <- unique(seurat_obj@meta.data[[col]])
    cat("\n  Unique values in '", col, "' (", length(vals), "):\n", sep = "")
    cat("   ", paste(head(vals, 30), collapse = ", "), "\n")
  } else if (!is.null(col) && col == "TODO") {
    cat("\n  WARNING: Column '", col, "' is still set to TODO — please update the config.\n", sep = "")
  }
}


# =============================================================================
# SECTION C: PSEUDOBULK AGGREGATION
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION C: Pseudobulk aggregation\n")
cat(strrep("=", 70), "\n\n")

library(Seurat)

# Validate columns exist
for (col in c(SAMPLE_COL, CELLTYPE_COL)) {
  if (col == "TODO") stop("Please set ", col, " in the config block at the top.")
  if (!col %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", col, "' not found in meta.data. Available: ",
         paste(colnames(seurat_obj@meta.data), collapse = ", "))
  }
}

# ── Cell counts per group ─────────────────────────────────────────────────────
group_key <- paste(seurat_obj@meta.data[[SAMPLE_COL]],
                   seurat_obj@meta.data[[CELLTYPE_COL]],
                   sep = "__")
cell_counts <- table(group_key)
cat("  Cell counts per pseudobulk group (sample__celltype):\n")
print(sort(cell_counts))

# ── Aggregate expression ──────────────────────────────────────────────────────
cat("\n  Running AggregateExpression ...\n")
group_by_cols <- c(SAMPLE_COL, CELLTYPE_COL)

agg_list <- AggregateExpression(
  seurat_obj,
  group.by    = group_by_cols,
  assays      = "RNA",
  slot        = "counts",
  return.seurat = FALSE
)

pb_counts <- agg_list[["RNA"]]   # genes × pseudobulk samples
cat("  Pseudobulk matrix dimensions (pre-filter):", nrow(pb_counts), "genes x",
    ncol(pb_counts), "samples\n")
cat("  Column names (first 10):", paste(head(colnames(pb_counts), 10), collapse = ", "), "\n")

# ── Parse column names into metadata ─────────────────────────────────────────
# AggregateExpression names columns as "sample_celltype" using the group.by columns
# The separator is "_" by default; we split on the last occurrence of CELLTYPE values
# to be robust, we use the known unique values
cat("\n  Parsing column names into sample/cell_type metadata...\n")

cell_types_known <- unique(seurat_obj@meta.data[[CELLTYPE_COL]])
samples_known    <- unique(seurat_obj@meta.data[[SAMPLE_COL]])

parse_pb_colname <- function(col_name, cell_types) {
  # Try to match a known cell type at the end (after the last "_")
  for (ct in sort(cell_types, decreasing = TRUE)) {
    # Seurat uses "_" as separator, cell type may contain spaces replaced by _
    ct_pattern <- gsub(" ", "_", ct)
    suffix     <- paste0("_", ct_pattern)
    if (endsWith(col_name, suffix)) {
      sample <- sub(paste0(suffix, "$"), "", col_name)
      return(list(sample = sample, cell_type = ct))
    }
  }
  # Fallback: split on last "_"
  parts <- strsplit(col_name, "_")[[1]]
  n <- length(parts)
  list(sample = paste(parts[seq_len(n - 1)], collapse = "_"),
       cell_type = parts[n])
}

pb_meta <- do.call(rbind, lapply(colnames(pb_counts), function(cn) {
  parsed <- parse_pb_colname(cn, cell_types_known)
  data.frame(pb_sample  = cn,
             sample     = parsed$sample,
             cell_type  = parsed$cell_type,
             stringsAsFactors = FALSE)
}))
rownames(pb_meta) <- pb_meta$pb_sample

if (!is.null(CONDITION_COL) && CONDITION_COL != "TODO" &&
    CONDITION_COL %in% colnames(seurat_obj@meta.data)) {
  # Map condition from sample -> condition via original metadata
  sample_cond_map <- seurat_obj@meta.data %>%
    dplyr::select(dplyr::all_of(c(SAMPLE_COL, CONDITION_COL))) %>%
    dplyr::distinct()
  rownames(sample_cond_map) <- sample_cond_map[[SAMPLE_COL]]
  pb_meta$condition <- sample_cond_map[pb_meta$sample, CONDITION_COL]
} else {
  pb_meta$condition <- NA
}

cat("  Pseudobulk metadata (first 5 rows):\n")
print(head(pb_meta, 5))

# ── Filter pseudobulk samples ─────────────────────────────────────────────────
cell_count_vec <- cell_counts[pb_meta$pb_sample]
total_counts_vec <- colSums(pb_counts)

keep_mask <- (cell_count_vec >= MIN_CELLS) & (total_counts_vec >= MIN_COUNTS)
cat("\n  Pre-filter samples:", ncol(pb_counts), "\n")
cat("  Removed (< MIN_CELLS =", MIN_CELLS, "cells):",
    sum(cell_count_vec < MIN_CELLS), "\n")
cat("  Removed (< MIN_COUNTS =", MIN_COUNTS, "counts):",
    sum(total_counts_vec < MIN_COUNTS), "\n")

pb_counts <- pb_counts[, keep_mask, drop = FALSE]
pb_meta   <- pb_meta[keep_mask, , drop = FALSE]
cat("  Post-filter samples:", ncol(pb_counts), "\n")


# =============================================================================
# SECTION D: NORMALIZE WITH EDGER
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION D: edgeR normalization (TMM + log-CPM)\n")
cat(strrep("=", 70), "\n\n")

dge <- DGEList(counts = pb_counts)
cat("  DGEList created:", nrow(dge), "genes x", ncol(dge), "samples\n")

# Filter lowly expressed genes
keep_genes <- filterByExpr(dge, min.count = 10)
cat("  Genes passing filterByExpr(min.count=10):", sum(keep_genes),
    "of", length(keep_genes), "\n")
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge, method = "TMM")
cat("  TMM normalization factors computed.\n")

# Log-CPM
log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  log-CPM matrix dimensions:", nrow(log_cpm), "genes x", ncol(log_cpm), "samples\n")


# =============================================================================
# SECTION E: DEFINE GENE SETS
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION E: Gene sets and coverage check\n")
cat(strrep("=", 70), "\n\n")

gene_sets <- list(
  NfkB = c("IKBKB", "NFKB1", "GFAP", "AQP4", "C3", "RELA", "RELB",
            "MAPK1", "MAPK3", "MAPK8", "IL6", "C1QA", "C1QB", "C1QC",
            "ICAM1", "TRAF1", "TRAF2", "CD44", "BCL2", "NOS2", "PTGS2",
            "LIF", "IL1B", "TNF", "CSF2", "CCL2", "CCL5", "CXCL1"),
  JAK_STAT = c("S100A10", "JAK1", "JAK2", "STAT1", "STAT3", "LIFR", "S100B",
               "SOCS3", "C1S", "C1R", "C3", "SERPINA3", "IL6", "CNTF",
               "IFNAR1", "IFNAR2"),
  WNT_BCATENIN = c("WNT1", "WNT3A", "WNT5A", "FZD1", "FZD2", "FZD9", "GSK3B",
                   "SOX9", "AXIN1", "SFRP1", "SFRP2", "FRZB", "TCF7", "LEF1",
                   "CD44", "SYT1", "SYT4", "CASP1", "CASP3", "CASP9", "GJA1", "GLUL"),
  PI3K_AKT = c("PIK3R1", "PIK3R2", "PIK3CA", "AKT1", "AKT2", "PTEN", "MTOR",
               "EGFR", "PDGFRA", "IGF1", "TLR2", "TLR4", "BDNF", "CDK1",
               "NFKB1", "BAD")
)

# ── Coverage check ────────────────────────────────────────────────────────────
all_expr_genes <- rownames(log_cpm)
coverage_rows  <- list()

for (pathway_name in names(gene_sets)) {
  gs         <- gene_sets[[pathway_name]]
  found      <- gs[gs %in% all_expr_genes]
  missing    <- gs[!gs %in% all_expr_genes]
  pct        <- 100 * length(found) / length(gs)

  cat("  Pathway:", pathway_name, "\n")
  cat("    Total genes:", length(gs), "\n")
  cat("    Found      :", length(found), "—", paste(found, collapse = ", "), "\n")
  cat("    Missing    :", length(missing), "—", paste(missing, collapse = ", "), "\n")
  cat("    Coverage   :", round(pct, 1), "%\n")
  if (pct < 50) {
    cat("    WARNING: Coverage < 50%! ssGSEA results for this pathway may be unreliable.\n")
  }
  cat("\n")

  coverage_rows[[pathway_name]] <- data.frame(
    pathway          = pathway_name,
    total_genes      = length(gs),
    found_genes      = length(found),
    missing_genes    = length(missing),
    coverage_pct     = round(pct, 1),
    found_list       = paste(found, collapse = ";"),
    missing_list     = paste(missing, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

coverage_df <- do.call(rbind, coverage_rows)
write.csv(coverage_df, "gsva_results/gene_coverage.csv", row.names = FALSE)
cat("  Coverage report saved to gsva_results/gene_coverage.csv\n")


# =============================================================================
# SECTION F: RUN GSVA / ssGSEA
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION F: Running ssGSEA\n")
cat(strrep("=", 70), "\n\n")

ssgsea_scores <- NULL

tryCatch({
  # GSVA >= 2.0 syntax
  cat("  Trying GSVA >= 2.0 syntax (ssgseaParam)...\n")
  params <- ssgseaParam(
    exprData  = log_cpm,
    geneSets  = gene_sets,
    normalize = TRUE
  )
  ssgsea_scores <- gsva(params,
                        BPPARAM = BiocParallel::SnowParam(workers = NCORES,
                                                          progressbar = TRUE))
  cat("  ssGSEA completed using GSVA >= 2.0 syntax.\n")
}, error = function(e) {
  cat("  GSVA >= 2.0 failed:", conditionMessage(e), "\n")
  cat("  Falling back to GSVA < 2.0 syntax...\n")
  tryCatch({
    ssgsea_scores <<- gsva(
      expr            = log_cpm,
      gset.idx.list   = gene_sets,
      method          = "ssgsea",
      ssgsea.norm     = TRUE,
      verbose         = TRUE,
      BPPARAM         = BiocParallel::SnowParam(workers = NCORES)
    )
    cat("  ssGSEA completed using GSVA < 2.0 syntax.\n")
  }, error = function(e2) {
    stop("Both GSVA syntaxes failed: ", conditionMessage(e2))
  })
})

cat("\n  ssGSEA score matrix dimensions:", nrow(ssgsea_scores), "pathways x",
    ncol(ssgsea_scores), "pseudobulk samples\n")
cat("  Score matrix head:\n")
print(head(ssgsea_scores))

write.csv(ssgsea_scores, "gsva_results/ssgsea_scores.csv")
cat("  Raw scores saved to gsva_results/ssgsea_scores.csv\n")


# =============================================================================
# SECTION G: VISUALIZATIONS
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION G: Generating visualizations\n")
cat(strrep("=", 70), "\n\n")

# ── Build annotation for heatmaps ──────────────────────────────────────────────
ann_col <- data.frame(
  cell_type = pb_meta$cell_type,
  sample    = pb_meta$sample,
  row.names = pb_meta$pb_sample
)
if (!all(is.na(pb_meta$condition))) {
  ann_col$condition <- pb_meta$condition
}
# Keep only columns present (and non-NA)
ann_col <- ann_col[, colSums(!is.na(ann_col)) > 0, drop = FALSE]

rdylbu_colors <- rev(brewer.pal(11, "RdYlBu"))

# ── G1: Heatmap (scaled rows) ─────────────────────────────────────────────────
cat("  G1: Heatmap (scaled) — heatmap_all_pathways.pdf\n")
pdf("gsva_results/heatmap_all_pathways.pdf", width = 14, height = 6)

pheatmap(
  ssgsea_scores,
  scale             = "row",
  annotation_col    = ann_col,
  clustering_method = "ward.D2",
  color             = rdylbu_colors,
  main              = "ssGSEA Scores (row-scaled)",
  fontsize          = 9,
  show_colnames     = ncol(ssgsea_scores) <= 50
)

# Second heatmap: raw (unscaled) scores
pheatmap(
  ssgsea_scores,
  scale             = "none",
  annotation_col    = ann_col,
  clustering_method = "ward.D2",
  color             = rdylbu_colors,
  main              = "ssGSEA Scores (raw, unscaled)",
  fontsize          = 9,
  show_colnames     = ncol(ssgsea_scores) <= 50
)
dev.off()
cat("  Saved.\n")

# ── G2: Boxplot by cell type ──────────────────────────────────────────────────
cat("  G2: Boxplot by cell type — boxplot_by_celltype.pdf\n")

score_long <- melt(as.data.frame(t(ssgsea_scores)),
                   variable.name = "pathway", value.name = "score")
score_long$pb_sample  <- rep(rownames(t(ssgsea_scores)), ncol(ssgsea_scores))
score_long$cell_type  <- pb_meta[score_long$pb_sample, "cell_type"]
score_long$sample     <- pb_meta[score_long$pb_sample, "sample"]
score_long$condition  <- pb_meta[score_long$pb_sample, "condition"]

p_box_ct <- ggplot(score_long, aes(x = pathway, y = score, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8)) +
  facet_wrap(~ pathway, scales = "free_y", nrow = 1) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "ssGSEA Scores by Cell Type", x = NULL, y = "ssGSEA Score",
       fill = "Cell Type")

pdf("gsva_results/boxplot_by_celltype.pdf", width = 16, height = 6)
print(p_box_ct)
dev.off()
cat("  Saved.\n")

# ── G3: Boxplot by condition (if available) ───────────────────────────────────
if (!is.null(CONDITION_COL) && CONDITION_COL != "TODO" &&
    !all(is.na(score_long$condition))) {
  cat("  G3: Boxplot by condition — boxplot_by_condition.pdf\n")

  p_box_cond <- ggplot(score_long,
                       aes(x = condition, y = score, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(aes(color = condition),
               position = position_jitterdodge(jitter.width = 0.2),
               size = 1.5, alpha = 0.7) +
    facet_wrap(~ pathway, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = "ssGSEA Scores by Condition", x = NULL, y = "ssGSEA Score",
         fill = "Condition", color = "Condition")

  pdf("gsva_results/boxplot_by_condition.pdf", width = 14, height = 8)
  print(p_box_cond)
  dev.off()
  cat("  Saved.\n")
} else {
  cat("  G3: Skipping — CONDITION_COL not set or no condition data available.\n")
}

# ── G4: Per-pathway dot plots ─────────────────────────────────────────────────
cat("  G4: Per-pathway dot plots — dotplot_pathways.pdf\n")

plot_list <- lapply(levels(score_long$pathway), function(pw) {
  df_pw <- score_long[score_long$pathway == pw, ]
  p <- ggplot(df_pw, aes(x = cell_type, y = score)) +
    geom_boxplot(aes(fill = cell_type), alpha = 0.3, outlier.shape = NA) +
    geom_point(aes(color = condition), size = 2, alpha = 0.8,
               position = position_jitter(width = 0.15)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(title = pw, x = NULL, y = "ssGSEA Score")
  p
})

pdf("gsva_results/dotplot_pathways.pdf", width = 12, height = 5)
for (p in plot_list) print(p)
dev.off()
cat("  Saved.\n")

# ── G5: Pathway-pathway scatter plots ─────────────────────────────────────────
cat("  G5: Scatter plots — scatter_pathway_pairs.pdf\n")

score_wide <- as.data.frame(t(ssgsea_scores))
score_wide$cell_type <- pb_meta$cell_type
score_wide$condition <- pb_meta$condition

pairs_to_plot <- list(
  c("NfkB",         "JAK_STAT"),
  c("NfkB",         "PI3K_AKT"),
  c("WNT_BCATENIN", "PI3K_AKT")
)

scatter_plots <- lapply(pairs_to_plot, function(pair) {
  pw1 <- pair[1]; pw2 <- pair[2]
  if (!all(c(pw1, pw2) %in% colnames(score_wide))) {
    cat("  Skipping scatter:", pw1, "vs", pw2, "— pathway not in scores.\n")
    return(NULL)
  }
  df_pair <- score_wide[, c(pw1, pw2, "cell_type", "condition")]
  colnames(df_pair)[1:2] <- c("x_score", "y_score")

  ct <- cor.test(df_pair$x_score, df_pair$y_score, method = "pearson")
  label <- sprintf("r = %.2f\np = %.3f", ct$estimate, ct$p.value)

  ggplot(df_pair, aes(x = x_score, y = y_score, color = cell_type)) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, label = label,
             hjust = -0.1, vjust = 1.2, size = 3.5) +
    theme_bw(base_size = 11) +
    labs(title = paste(pw1, "vs", pw2),
         x = paste(pw1, "ssGSEA Score"),
         y = paste(pw2, "ssGSEA Score"),
         color = "Cell Type")
})

scatter_plots <- Filter(Negate(is.null), scatter_plots)

pdf("gsva_results/scatter_pathway_pairs.pdf", width = 5, height = 4.5)
for (p in scatter_plots) print(p)
dev.off()
cat("  Saved.\n")

# ── G6: Pathway-pathway correlation heatmap ───────────────────────────────────
cat("  G6: Pathway correlation heatmap — pathway_correlation.pdf\n")

pathway_cor <- cor(t(ssgsea_scores), method = "pearson")

pdf("gsva_results/pathway_correlation.pdf", width = 7, height = 6)
pheatmap(
  pathway_cor,
  display_numbers = TRUE,
  number_format   = "%.2f",
  color           = rdylbu_colors,
  main            = "Pathway-Pathway Pearson Correlation",
  fontsize        = 11,
  breaks          = seq(-1, 1, length.out = length(rdylbu_colors) + 1)
)
dev.off()
cat("  Saved.\n")


# =============================================================================
# SECTION H: STATISTICAL TESTING WITH LIMMA
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION H: Differential enrichment testing with limma\n")
cat(strrep("=", 70), "\n\n")

all_limma_results <- list()

if (is.null(CONDITION_COL) || CONDITION_COL == "TODO" ||
    all(is.na(pb_meta$condition))) {
  cat("  CONDITION_COL not set — skipping limma testing.\n")
} else {

  cell_types_pb <- unique(pb_meta$cell_type)
  cat("  Cell types to test:", paste(cell_types_pb, collapse = ", "), "\n\n")

  for (ct in cell_types_pb) {
    cat("  Testing cell type:", ct, "\n")

    ct_mask  <- pb_meta$cell_type == ct
    ct_scores <- ssgsea_scores[, ct_mask, drop = FALSE]
    ct_meta   <- pb_meta[ct_mask, ]

    if (ncol(ct_scores) < 3) {
      cat("  Skipping — fewer than 3 pseudobulk samples.\n\n")
      next
    }

    cond_vec <- ct_meta$condition
    cond_tab <- table(cond_vec)
    if (length(cond_tab) < 2) {
      cat("  Skipping — only one condition level:", names(cond_tab), "\n\n")
      next
    }
    if (any(cond_tab < 2)) {
      cat("  Skipping — a condition has fewer than 2 samples.\n\n")
      next
    }

    cat("  Condition counts:", paste(names(cond_tab), cond_tab, sep = "=",
                                      collapse = ", "), "\n")

    tryCatch({
      cond_factor <- factor(cond_vec)
      design      <- model.matrix(~ cond_factor)
      colnames(design) <- make.names(colnames(design))

      fit   <- lmFit(ct_scores, design)
      fit   <- eBayes(fit)

      # Test the condition coefficient (second column)
      coef_name <- colnames(design)[2]
      tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
      tt$pathway   <- rownames(tt)
      tt$cell_type <- ct
      all_limma_results[[ct]] <- tt

      cat("  Significant pathways (adj.P.Val < 0.05):",
          sum(tt$adj.P.Val < 0.05, na.rm = TRUE), "\n")
    }, error = function(e) {
      cat("  ERROR during limma for", ct, ":", conditionMessage(e), "\n")
    })
    cat("\n")
  }
}

if (length(all_limma_results) > 0) {
  combined_results <- do.call(rbind, all_limma_results)
  rownames(combined_results) <- NULL
  write.csv(combined_results, "gsva_results/limma_differential_enrichment.csv",
            row.names = FALSE)
  cat("  limma results saved to gsva_results/limma_differential_enrichment.csv\n")

  sig_results <- combined_results[!is.na(combined_results$adj.P.Val) &
                                    combined_results$adj.P.Val < 0.05, ]
  cat("\n  Significant results (adj.P.Val < 0.05):\n")
  if (nrow(sig_results) == 0) {
    cat("  None found.\n")
  } else {
    print(sig_results[, c("pathway", "cell_type", "logFC", "adj.P.Val")])
  }
} else {
  cat("  No limma results generated.\n")
  combined_results <- data.frame()
}


# =============================================================================
# SECTION I: FINAL SUMMARY AND SAVE OUTPUTS
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  SECTION I: Final summary and saving outputs\n")
cat(strrep("=", 70), "\n\n")

cat("  ── SUMMARY ──────────────────────────────────────────────\n")
cat("  Total pseudobulk samples (post-filter):", ncol(pb_counts), "\n")
cat("  Number of cell types:", length(unique(pb_meta$cell_type)), "\n")
cat("  Cell types:", paste(unique(pb_meta$cell_type), collapse = ", "), "\n")
cat("  Genes passing edgeR filter:", nrow(dge), "\n")
cat("  Pathways analyzed:", nrow(ssgsea_scores), "\n")

cat("\n  Gene coverage per pathway:\n")
for (i in seq_len(nrow(coverage_df))) {
  cat(sprintf("    %-20s : %d/%d (%.1f%%)\n",
              coverage_df$pathway[i], coverage_df$found_genes[i],
              coverage_df$total_genes[i], coverage_df$coverage_pct[i]))
}

if (nrow(combined_results) > 0) {
  sig_n <- sum(combined_results$adj.P.Val < 0.05, na.rm = TRUE)
  cat("\n  Significant differential enrichment results (adj.P.Val < 0.05):", sig_n, "\n")
} else {
  cat("\n  No differential enrichment testing performed.\n")
}

# Save DGEList for future use
saveRDS(dge, "gsva_results/dge_pseudobulk.rds")
cat("\n  DGEList saved to gsva_results/dge_pseudobulk.rds\n")

# ── List all output files ─────────────────────────────────────────────────────
cat("\n  ── OUTPUT FILES ─────────────────────────────────────────\n")
output_files <- list.files("gsva_results", full.names = TRUE)
for (f in sort(output_files)) {
  size_kb <- round(file.info(f)$size / 1024, 1)
  cat(sprintf("    %-55s  %6.1f KB\n", f, size_kb))
}

cat("\n", strrep("=", 70), "\n")
cat("  Pipeline complete. Results in gsva_results/\n")
cat(strrep("=", 70), "\n\n")
