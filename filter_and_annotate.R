#!/usr/bin/env Rscript
# =============================================================================
# filter_and_annotate.R
# =============================================================================
# Filter cells and annotate cell types for PPMI iPSC scRNA-seq data.
#
# Author : Manveer Chauhan
# Usage  : Rscript filter_and_annotate.R
#
# Input  : export/seurat_object.rds
# Output : export/seurat_filtered.rds
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_RDS  <- "export/seurat_object.rds"
OUTPUT_RDS <- "export/seurat_filtered.rds"

DONOR_COL  <- "PPMI Donor ID"   # column used to remove doublets/unassigned

# =============================================================================
# HELPER
# =============================================================================
section <- function(title) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  ", title, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

cell_counts <- function(obj, label = "") {
  cat(sprintf("  %-35s %d cells\n", label, ncol(obj)))
}

# =============================================================================
# LOAD
# =============================================================================
section("LOADING SEURAT OBJECT")

library(Seurat)
cat("  Seurat version:", as.character(packageVersion("Seurat")), "\n")

if (!file.exists(INPUT_RDS)) {
  stop("Could not find: ", INPUT_RDS, "\n  Run h5ad_to_seurat.R first.")
}
seurat_obj <- readRDS(INPUT_RDS)
cat("  Loaded:", INPUT_RDS, "\n")
n_cells_original <- ncol(seurat_obj)
cell_counts(seurat_obj, "Cells at load:")

# =============================================================================
# SECTION 1: REMOVE DOUBLETS AND UNASSIGNED CELLS
# =============================================================================
section("SECTION 1: REMOVE DOUBLETS & UNASSIGNED CELLS")

donor_vals <- seurat_obj@meta.data[[DONOR_COL]]

# Summarise before filtering
cat("\n  PPMI Donor ID value counts (pre-filter):\n")
vc <- sort(table(donor_vals), decreasing = TRUE)
for (i in seq_along(vc)) {
  cat(sprintf("    %-20s  %d\n", names(vc)[i], vc[i]))
}

# Identify cells to remove
to_remove <- donor_vals %in% c("doublet", "unassigned")
n_doublet    <- sum(donor_vals == "doublet",    na.rm = TRUE)
n_unassigned <- sum(donor_vals == "unassigned", na.rm = TRUE)
cat(sprintf("\n  Removing %d doublets + %d unassigned = %d cells total\n",
            n_doublet, n_unassigned, sum(to_remove)))

seurat_obj <- seurat_obj[, !to_remove]
cell_counts(seurat_obj, "Cells after doublet/unassigned removal:")

# =============================================================================
# SECTION 2: ADD CONDITION METADATA
# =============================================================================
section("SECTION 2: ASSIGN CONDITION LABELS")

donor_to_condition <- c(
  PPMISI4106  = "LRRK2_PD",  PPMISI50860 = "LRRK2_PD",
  PPMISI40273 = "LRRK2_PD",  PPMISI51440 = "LRRK2_PD",
  PPMISI51625 = "LRRK2_PD",  PPMISI51330 = "LRRK2_PD",
  PPMISI3220  = "Idiopathic_PD", PPMISI4102  = "Idiopathic_PD",
  PPMISI3473  = "Idiopathic_PD", PPMISI4099  = "Idiopathic_PD",
  PPMISI3448  = "Idiopathic_PD",
  PPMISI3952  = "Control", PPMISI3452  = "Control",
  PPMISI3966  = "Control", PPMISI4105  = "Control",
  PPMISI3411  = "Control"
)

# Check for any donor IDs in the object not covered by the mapping
donors_in_obj <- unique(seurat_obj@meta.data[[DONOR_COL]])
unmapped <- donors_in_obj[!donors_in_obj %in% names(donor_to_condition)]
if (length(unmapped) > 0) {
  stop("Unmapped donor ID(s) found — update donor_to_condition:\n  ",
       paste(unmapped, collapse = ", "))
}

# Map and add to metadata
seurat_obj@meta.data$condition <- unname(donor_to_condition[seurat_obj@meta.data[[DONOR_COL]]])
cat("  'condition' column added.\n\n")

# Verification table: cells per condition x donor
cat("  Cell counts per condition:\n")
tbl <- table(Condition = seurat_obj@meta.data$condition)
print(tbl)

cat("\n  Cell counts per donor x condition:\n")
tbl2 <- table(Donor     = seurat_obj@meta.data[[DONOR_COL]],
              Condition = seurat_obj@meta.data$condition)
print(tbl2)

# =============================================================================
# SECTION 3: REMOVE LRRK2_PD CELLS
# =============================================================================
section("SECTION 3: REMOVE LRRK2_PD CELLS")

n_before_lrrk2 <- ncol(seurat_obj)
seurat_obj <- seurat_obj[, seurat_obj@meta.data$condition != "LRRK2_PD"]
n_removed_lrrk2 <- n_before_lrrk2 - ncol(seurat_obj)
cat(sprintf("  Removed %d LRRK2_PD cells.\n", n_removed_lrrk2))
cell_counts(seurat_obj, "Cells retained (Idiopathic_PD + Control):")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
section("SUMMARY")

cat("\n  Object:\n")
print(seurat_obj)

cat("\n  Metadata columns:\n")
for (col in colnames(seurat_obj@meta.data)) {
  n_unique <- length(unique(seurat_obj@meta.data[[col]]))
  cat(sprintf("    %-30s  n_unique = %d\n", col, n_unique))
}

# =============================================================================
# SAVE
# =============================================================================
section("SAVING")

cat("  Saving to:", OUTPUT_RDS, "...\n")
saveRDS(seurat_obj, file = OUTPUT_RDS)
cat("  Saved. File size:", round(file.size(OUTPUT_RDS) / 1e6, 1), "MB\n")

cat(sprintf("
  INPUT  : %s
  OUTPUT : %s
  Cells retained : %d / %d

  To load:
    seurat_obj <- readRDS(\"%s\")

", INPUT_RDS, OUTPUT_RDS,
   ncol(seurat_obj),
   n_cells_original,
   OUTPUT_RDS))
