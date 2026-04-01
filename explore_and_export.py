#!/usr/bin/env python3
"""
explore_and_export.py
=====================
Part 1: Explore a scRNA-seq AnnData (.h5ad) object and export data into
R-friendly formats for downstream pseudobulk ssGSEA analysis.

Author: Manveer Chauhan
Usage:  python explore_and_export.py

Outputs (written to export/):
  - adata_rawcounts.h5ad   : AnnData with raw counts in .X
  - matrix.mtx             : Sparse count matrix (genes x cells)
  - barcodes.tsv           : Cell barcodes
  - genes.tsv              : Gene names
  - metadata.csv           : Full .obs dataframe
"""

import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
from pathlib import Path

# ── Try to import scanpy / anndata ────────────────────────────────────────────
try:
    import scanpy as sc
    import anndata as ad
except ImportError:
    raise ImportError(
        "scanpy and anndata are required. Install with:\n"
        "  pip install scanpy anndata"
    )

# =============================================================================
# CONFIGURATION
# =============================================================================
H5AD_FILE   = "PPMI_156_MJFF-1651.scRNAseq.data.h5ad"
EXPORT_DIR  = Path("export")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def section(title: str) -> None:
    """Print a clearly visible section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)


def check_integer_matrix(X) -> dict:
    """
    Return dtype, min, max, and whether all values are integers.
    Works for dense arrays and scipy sparse matrices.
    """
    if scipy.sparse.issparse(X):
        data = X.data  # only non-zero values
        if len(data) == 0:
            return {"dtype": X.dtype, "min": 0, "max": 0, "is_integer": True,
                    "sparsity": 1.0}
        vmin = float(data.min())
        vmax = float(data.max())
        is_int = bool(np.all(data == np.floor(data)))
        sparsity = 1.0 - (X.nnz / (X.shape[0] * X.shape[1]))
    else:
        X_arr = np.asarray(X)
        vmin = float(X_arr.min())
        vmax = float(X_arr.max())
        is_int = bool(np.all(X_arr == np.floor(X_arr)))
        sparsity = float(np.mean(X_arr == 0))
    return {
        "dtype": str(X.dtype),
        "min": vmin,
        "max": vmax,
        "is_integer": is_int,
        "sparsity": sparsity,
    }


# =============================================================================
# LOAD DATA
# =============================================================================
section("LOADING H5AD FILE")
print(f"  File: {H5AD_FILE}")
if not os.path.exists(H5AD_FILE):
    raise FileNotFoundError(
        f"Could not find {H5AD_FILE}\n"
        "  Please place the file in the current working directory."
    )
adata = sc.read_h5ad(H5AD_FILE)
print(f"  Loaded successfully.")


# =============================================================================
# SECTION 1: BASIC INFO
# =============================================================================
section("SECTION 1: BASIC OBJECT INFO")

print(f"\n  Shape (cells × genes): {adata.shape}")
print(f"  n_obs (cells)        : {adata.n_obs:,}")
print(f"  n_vars (genes)       : {adata.n_vars:,}")

print(f"\n  .obs columns ({len(adata.obs.columns)}):")
for col in adata.obs.columns.tolist():
    print(f"    {col}")

print(f"\n  .var head(10):")
print(adata.var.head(10).to_string())

print(f"\n  .var_names[:20]:")
print(list(adata.var_names[:20]))


# =============================================================================
# SECTION 2: .obs COLUMN INSPECTION
# =============================================================================
section("SECTION 2: .obs COLUMN INSPECTION")

for col in adata.obs.columns:
    dtype   = adata.obs[col].dtype
    n_unique = adata.obs[col].nunique()
    uniq_vals = list(adata.obs[col].unique()[:20])

    print(f"\n  Column : {col}")
    print(f"    dtype    : {dtype}")
    print(f"    n_unique : {n_unique}")
    print(f"    first 20 unique values: {uniq_vals}")

    if hasattr(dtype, "name") and (
        dtype.name == "category" or str(dtype) == "object"
    ):
        vc = adata.obs[col].value_counts().head(20)
        print(f"    value_counts (top 20):")
        print(vc.to_string(header=False))


# =============================================================================
# SECTION 3: LAYERS INSPECTION
# =============================================================================
section("SECTION 3: LAYERS INSPECTION")

if len(adata.layers) == 0:
    print("  No layers found.")
else:
    print(f"  Found {len(adata.layers)} layer(s): {list(adata.layers.keys())}")
    for layer_name, layer_data in adata.layers.items():
        info = check_integer_matrix(layer_data)
        print(f"\n  Layer: '{layer_name}'")
        print(f"    dtype      : {info['dtype']}")
        print(f"    min        : {info['min']}")
        print(f"    max        : {info['max']}")
        print(f"    is_integer : {info['is_integer']}")
        print(f"    sparsity   : {info['sparsity']:.3f}")
        if info["is_integer"] and info["min"] >= 0:
            print(f"    >>> Looks like RAW COUNTS")
        else:
            print(f"    >>> Likely NORMALIZED data")


# =============================================================================
# SECTION 4: .X INSPECTION
# =============================================================================
section("SECTION 4: .X INSPECTION")

x_info = check_integer_matrix(adata.X)
print(f"  dtype      : {x_info['dtype']}")
print(f"  min        : {x_info['min']}")
print(f"  max        : {x_info['max']}")
print(f"  is_integer : {x_info['is_integer']}")
print(f"  sparsity   : {x_info['sparsity']:.3f}")

if x_info["is_integer"] and x_info["min"] >= 0:
    x_verdict = "RAW COUNTS"
else:
    x_verdict = "NORMALIZED DATA (not raw counts)"
print(f"\n  .X verdict: {x_verdict}")


# =============================================================================
# SECTION 5: .raw INSPECTION
# =============================================================================
section("SECTION 5: .raw INSPECTION")

has_raw = adata.raw is not None
print(f"  adata.raw exists: {has_raw}")

raw_is_integer = False
if has_raw:
    print(f"\n  adata.raw.var_names[:20]:")
    print(list(adata.raw.var_names[:20]))
    raw_info = check_integer_matrix(adata.raw.X)
    print(f"\n  adata.raw.X:")
    print(f"    dtype      : {raw_info['dtype']}")
    print(f"    min        : {raw_info['min']}")
    print(f"    max        : {raw_info['max']}")
    print(f"    is_integer : {raw_info['is_integer']}")
    print(f"    sparsity   : {raw_info['sparsity']:.3f}")
    raw_is_integer = raw_info["is_integer"] and raw_info["min"] >= 0
    if raw_is_integer:
        print("  >>> .raw.X looks like RAW COUNTS")
    else:
        print("  >>> .raw.X does NOT look like raw counts")


# =============================================================================
# SECTION 6: .obsm EMBEDDINGS
# =============================================================================
section("SECTION 6: .obsm EMBEDDINGS")

if len(adata.obsm) == 0:
    print("  No embeddings found in .obsm.")
else:
    print(f"  Found {len(adata.obsm)} embedding(s):")
    for key in adata.obsm.keys():
        shape = adata.obsm[key].shape
        print(f"    {key} — shape: {shape}")


# =============================================================================
# SECTION 7: RECOMMENDATION
# =============================================================================
section("SECTION 7: RECOMMENDATION — WHERE ARE THE RAW COUNTS?")

# Determine the best source of raw counts
raw_count_source = None
raw_count_layer  = None

# Priority 1: .X is already raw counts
if x_info["is_integer"] and x_info["min"] >= 0:
    raw_count_source = ".X"
    print("  RECOMMENDATION: Raw counts appear to be in .X")

# Priority 2: .raw.X has integer counts (and .X is normalized)
elif has_raw and raw_is_integer:
    raw_count_source = ".raw.X"
    print("  RECOMMENDATION: Raw counts appear to be in .raw.X")
    print("  ACTION: Will convert .raw to AnnData before export")

# Priority 3: a layer holds integer counts
else:
    integer_layers = []
    for lname, ldata in adata.layers.items():
        info = check_integer_matrix(ldata)
        if info["is_integer"] and info["min"] >= 0:
            integer_layers.append(lname)

    if integer_layers:
        raw_count_layer  = integer_layers[0]
        raw_count_source = f"layers['{raw_count_layer}']"
        print(f"  RECOMMENDATION: Raw counts appear to be in layer '{raw_count_layer}'")
        if len(integer_layers) > 1:
            print(f"  (Other integer layers found: {integer_layers[1:]})")
        print(f"  ACTION: Will copy layer '{raw_count_layer}' to .X before export")
    else:
        raw_count_source = ".X (UNCERTAIN — no clear integer counts found)"
        print("  WARNING: Could not definitively identify raw counts.")
        print("  Defaulting to .X — check the output carefully.")

print(f"\n  >> Set COUNTS_LAYER in pseudobulk_gsva.R to: "
      f"{'NULL' if raw_count_layer is None else repr(raw_count_layer)}")


# =============================================================================
# SECTION 8: EXPORT
# =============================================================================
section("SECTION 8: EXPORTING DATA")

EXPORT_DIR.mkdir(parents=True, exist_ok=True)
print(f"  Export directory: {EXPORT_DIR.resolve()}")

# ── Build the export AnnData ──────────────────────────────────────────────────
print("\n  Preparing export AnnData...")

if raw_count_source == ".raw.X":
    # Convert .raw to a proper AnnData (genes × cells filtered to raw genes)
    adata_export = adata.raw.to_adata()
    print("  Using adata.raw.to_adata() as export object")

elif raw_count_layer is not None:
    # Copy the integer layer into .X
    import copy
    adata_export = adata.copy()
    adata_export.X = adata_export.layers[raw_count_layer].copy()
    # Optionally drop layers to save space
    adata_export.layers = {}
    print(f"  Copied layer '{raw_count_layer}' to .X")

else:
    # .X already holds the best counts (or we're uncertain — use .X as-is)
    adata_export = adata.copy()
    print("  Using .X as-is for export")

# ── Save adata_rawcounts.h5ad ─────────────────────────────────────────────────
h5ad_out = EXPORT_DIR / "adata_rawcounts.h5ad"
print(f"\n  Saving {h5ad_out} ...")
adata_export.write_h5ad(h5ad_out)
print(f"  Saved. Shape: {adata_export.shape}")

# ── Get count matrix for MTX export ──────────────────────────────────────────
X_export = adata_export.X
if not scipy.sparse.issparse(X_export):
    X_export = scipy.sparse.csc_matrix(X_export)
else:
    X_export = X_export.tocsc()

# scipy.io.mmwrite expects (features × barcodes) = (genes × cells)
# adata stores (cells × genes), so transpose for MTX convention
X_export_T = X_export.T   # genes × cells

# ── matrix.mtx ───────────────────────────────────────────────────────────────
mtx_out = EXPORT_DIR / "matrix.mtx"
print(f"\n  Saving {mtx_out} ...")
scipy.io.mmwrite(str(mtx_out), X_export_T)
print(f"  Saved. Shape (genes × cells): {X_export_T.shape}")

# ── barcodes.tsv ─────────────────────────────────────────────────────────────
barcodes_out = EXPORT_DIR / "barcodes.tsv"
print(f"\n  Saving {barcodes_out} ...")
pd.Series(adata_export.obs_names).to_csv(barcodes_out, index=False, header=False)
print(f"  Saved. {len(adata_export.obs_names):,} barcodes.")

# ── genes.tsv ────────────────────────────────────────────────────────────────
genes_out = EXPORT_DIR / "genes.tsv"
print(f"\n  Saving {genes_out} ...")
pd.Series(adata_export.var_names).to_csv(genes_out, index=False, header=False)
print(f"  Saved. {len(adata_export.var_names):,} genes.")

# ── metadata.csv ─────────────────────────────────────────────────────────────
meta_out = EXPORT_DIR / "metadata.csv"
print(f"\n  Saving {meta_out} ...")
# Use the original adata.obs (full metadata), aligned to export barcodes
obs_export = adata.obs.loc[adata_export.obs_names] if raw_count_source == ".raw.X" else adata_export.obs
obs_export.to_csv(meta_out)
print(f"  Saved. Shape: {obs_export.shape}")


# =============================================================================
# SECTION 9: FINAL SUMMARY
# =============================================================================
section("SECTION 9: FINAL SUMMARY")

print(f"""
  INPUT FILE   : {H5AD_FILE}
  OBJECT SHAPE : {adata.shape}  (cells × genes)

  RAW COUNT LOCATION : {raw_count_source}
  COUNTS_LAYER (R)   : {'NULL' if raw_count_layer is None else repr(raw_count_layer)}

  EXPORT FILES (in {EXPORT_DIR}/):
    adata_rawcounts.h5ad  — AnnData with raw counts in .X
    matrix.mtx            — Sparse matrix (genes × cells)
    barcodes.tsv          — {adata_export.n_obs:,} cell barcodes
    genes.tsv             — {adata_export.n_vars:,} genes
    metadata.csv          — {obs_export.shape[1]} metadata columns

  NEXT STEPS:
    1. Review SECTION 2 output above to identify:
         SAMPLE_COL    — column for biological replicates / donors
         CELLTYPE_COL  — column for cell type annotations
         CONDITION_COL — column for disease vs. control (or set NULL)
    2. Open pseudobulk_gsva.R and fill in the TODO config block.
    3. Run:  Rscript pseudobulk_gsva.R
""")

print("  Done.")
