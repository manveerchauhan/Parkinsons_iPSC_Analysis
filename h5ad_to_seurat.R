#!/usr/bin/env Rscript
# =============================================================================
# h5ad_to_seurat.R
# =============================================================================
# Convert an AnnData (.h5ad) file to a Seurat object and save as .rds.
#
# Author : Manveer Chauhan
# Usage  : Rscript h5ad_to_seurat.R
#
# Input  : export/adata_rawcounts.h5ad
# Output : export/seurat_object.rds
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================
MTX_PATH      <- "export/matrix.mtx"
BARCODES_PATH <- "export/barcodes.tsv"
GENES_PATH    <- "export/genes.tsv"
METADATA_PATH <- "export/metadata.csv"
OUTPUT_RDS    <- "export/seurat_object.rds"

# =============================================================================
# HELPER
# =============================================================================
section <- function(title) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  ", title, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

# =============================================================================
# INSTALL / LOAD PACKAGES
# =============================================================================
section("INSTALLING / LOADING PACKAGES")

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat", repos = "https://cloud.r-project.org")
}
library(Seurat)
cat("  Seurat loaded (version", as.character(packageVersion("Seurat")), ")\n")

# =============================================================================
# CHECK INPUT FILES
# =============================================================================
section("CHECKING INPUT FILES")

required_files <- c(MTX_PATH, BARCODES_PATH, GENES_PATH, METADATA_PATH)
for (f in required_files) {
  if (!file.exists(f)) {
    stop("Could not find: ", f, "\n  Run explore_and_export.py first.")
  }
  cat("  Found:", f, "\n")
}

# =============================================================================
# BUILD SEURAT OBJECT FROM MTX EXPORT
# =============================================================================
section("BUILDING SEURAT OBJECT FROM MTX FILES")

# ReadMtx expects: mtx (genes x cells), cells file, features file
cat("  Reading MTX files...\n")
counts <- ReadMtx(
  mtx      = MTX_PATH,
  cells    = BARCODES_PATH,
  features = GENES_PATH,
  feature.column = 1   # genes.tsv has one column (gene name)
)
cat(sprintf("  Count matrix: %d genes x %d cells\n", nrow(counts), ncol(counts)))

cat("  Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts  = counts,
  project = "PPMI_iPSC"
)
rm(counts); gc()
cat("  Seurat object created.\n")

# ── Attach metadata ───────────────────────────────────────────────────────────
cat("  Attaching metadata from", METADATA_PATH, "...\n")
meta <- read.csv(METADATA_PATH, row.names = 1, check.names = FALSE)

# Align to cells in the Seurat object
shared_cells <- intersect(colnames(seurat_obj), rownames(meta))
cat(sprintf("  Cells in Seurat: %d | Cells in metadata: %d | Shared: %d\n",
            ncol(seurat_obj), nrow(meta), length(shared_cells)))

for (col in colnames(meta)) {
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta[[col]][match(colnames(seurat_obj), rownames(meta))], col.name = col)
}
cat("  Metadata columns added:", paste(colnames(meta), collapse = ", "), "\n")

# =============================================================================
# INSPECT SEURAT OBJECT
# =============================================================================
section("SEURAT OBJECT SUMMARY")

cat("\n  Object:\n")
print(seurat_obj)

cat("\n  Assay(s)     :", paste(Assays(seurat_obj), collapse = ", "), "\n")
cat("  Active assay :", DefaultAssay(seurat_obj), "\n")
cat("  Cells        :", ncol(seurat_obj), "\n")
cat("  Features     :", nrow(seurat_obj), "\n")

cat("\n  Metadata columns:\n")
for (col in colnames(seurat_obj@meta.data)) {
  n_unique <- length(unique(seurat_obj@meta.data[[col]]))
  cat(sprintf("    %-30s  n_unique = %d\n", col, n_unique))
}

reductions <- Reductions(seurat_obj)
if (length(reductions) > 0) {
  cat("\n  Reductions   :", paste(reductions, collapse = ", "), "\n")
} else {
  cat("\n  Reductions   : none\n")
}

# =============================================================================
# SAVE
# =============================================================================
section("SAVING SEURAT OBJECT")

dir.create(dirname(OUTPUT_RDS), showWarnings = FALSE, recursive = TRUE)
cat("  Saving to:", OUTPUT_RDS, "...\n")
saveRDS(seurat_obj, file = OUTPUT_RDS)
cat("  Saved. File size:", round(file.size(OUTPUT_RDS) / 1e6, 1), "MB\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
section("DONE")

cat(sprintf("
  INPUT  : %s (+ barcodes, genes, metadata)
  OUTPUT : %s

  Cells    : %d
  Features : %d
  Metadata : %s

  To load in a future R session:
    seurat_obj <- readRDS(\"%s\")

", MTX_PATH, OUTPUT_RDS,
   ncol(seurat_obj), nrow(seurat_obj),
   paste(colnames(seurat_obj@meta.data), collapse = ", "),
   OUTPUT_RDS))
