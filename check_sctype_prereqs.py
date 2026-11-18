"""
check_sctype_prereqs.py
Run this FIRST before building the real scType runner.

scType (via sctypepy, a Python port) is marker-gene-based - no trained
classifier, no reference dataset, no tissue-mismatch risk like the
immune atlas had. It ships a built-in marker database that specifically
includes a "Liver" tissue category, which is directly relevant here.

This checks:
  1. Is sctypepy installed? What does its API actually look like?
  2. Does it run successfully on a small subset of the real data?
  3. What does the output actually contain (labels? scores? shape?) -
     needed to know how to extract a real confidence value, not just
     assume the documented API matches exactly.
"""

import sys

print("=" * 70)
print("1. sctypepy installation check")
print("=" * 70)
try:
    import sctypepy
    print(f"  sctypepy imported successfully")
    print(f"  Available functions: {[x for x in dir(sctypepy) if not x.startswith('_')]}")
except ImportError:
    print("  NOT INSTALLED. Run: pip install sctypepy")
    sys.exit(1)

print("\n" + "=" * 70)
print("2. Available tissue types in the built-in marker database")
print("=" * 70)
try:
    from sctypepy import get_available_tissues
    tissues = get_available_tissues()
    print(f"  {tissues}")
    if "Liver" not in str(tissues):
        print("  WARNING: 'Liver' not found in available tissues - check exact naming.")
except Exception as e:
    print(f"  Could not call get_available_tissues(): {e}")

print("\n" + "=" * 70)
print("3. Test run on a small subset of real data")
print("=" * 70)
import scanpy as sc
from sctypepy import run_sctype

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v6.h5ad"
print(f"Loading {H5AD_PATH} (subsetting to 2000 cells for a quick test)...")
adata = sc.read_h5ad(H5AD_PATH)
adata_test = adata[:2000].copy()

# scType typically expects log-normalized data in .X
adata_test.X = adata_test.layers["log1p_norm"].copy()

try:
    result = run_sctype(adata_test, tissue_type="Liver", groupby="leiden_r05")
    print("  run_sctype() succeeded")
    print(f"\n  obs columns added: {[c for c in result.obs.columns if 'sctype' in c.lower()]}")
    print(f"  obsm keys added: {[k for k in result.obsm.keys() if 'sctype' in k.lower()]}")

    if "sctype_classification" in result.obs.columns:
        print(f"\n  Sample classifications:\n{result.obs['sctype_classification'].value_counts()}")

    for key in result.obsm.keys():
        if "sctype" in key.lower():
            scores = result.obsm[key]
            print(f"\n  obsm['{key}'] shape: {scores.shape}, type: {type(scores)}")
            if hasattr(scores, "columns"):
                print(f"  columns: {list(scores.columns)[:10]}")
            print(f"  sample values (first row): {scores.iloc[0] if hasattr(scores, 'iloc') else scores[0]}")

except Exception as e:
    print(f"  ERROR running run_sctype(): {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("Paste this full output back - it tells us exactly how to extract")
print("labels + a real confidence value for the actual runner.")
print("=" * 70)