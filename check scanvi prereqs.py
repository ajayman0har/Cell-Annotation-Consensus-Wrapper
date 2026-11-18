"""
check_scanvi_prereqs.py
Run this FIRST, before we write any actual scANVI training code.

scVI/scANVI (via scvi-tools) is architecturally different from
CellTypist/SingleR: it's a variational autoencoder that needs to be
trained on RAW COUNTS (not log-normalized data - it does its own internal
normalization as part of the generative model), and training speed
depends heavily on whether you have a GPU available.

This checks:
  1. Is scvi-tools installed? What version?
  2. Is a GPU available (dramatically affects feasible training time)?
  3. Does the Guilliams reference have a raw counts layer?
  4. Does the Xiong query data have a raw counts layer?
  5. How many genes are shared between reference and query (scANVI needs
     a consistent gene set)?
"""

import sys

print("=" * 70)
print("1. scvi-tools installation check")
print("=" * 70)
try:
    import scvi
    print(f"  scvi-tools version: {scvi.__version__}")
except ImportError:
    print("  NOT INSTALLED. Run: pip install scvi-tools")
    sys.exit(1)

print("\n" + "=" * 70)
print("2. GPU availability check")
print("=" * 70)
import torch
gpu_available = torch.cuda.is_available()
print(f"  CUDA GPU available: {gpu_available}")
if gpu_available:
    print(f"  GPU: {torch.cuda.get_device_name(0)}")
else:
    print("  No GPU detected - training will run on CPU, which is "
          "meaningfully slower for scVI/scANVI than for CellTypist/SingleR. "
          "This affects how many epochs/how large a reference we can "
          "practically use.")

print("\n" + "=" * 70)
print("3. Guilliams reference - raw counts check")
print("=" * 70)
import scanpy as sc
GUILLIAMS_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\guilliams_output\guilliams_stst_50k_annotated.h5ad"
ref = sc.read_h5ad(GUILLIAMS_PATH)
print(f"  {ref.n_obs} cells x {ref.n_vars} genes")
print(f"  Layers available: {list(ref.layers.keys())}")
if "counts" in ref.layers:
    print(f"  'counts' layer max value: {ref.layers['counts'].max()} "
          f"(should look like raw integer counts, e.g. tens to hundreds)")
else:
    print("  WARNING: no 'counts' layer found - scANVI needs raw counts, "
          "not the log1p-normalised .X")
print(f"  Label column: 'cell_type', {ref.obs['cell_type'].nunique()} categories")

print("\n" + "=" * 70)
print("4. Xiong query data - raw counts check")
print("=" * 70)
QUERY_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\Scripts\xiong_immune_output\xiong_wrapper_v2.h5ad"
query = sc.read_h5ad(QUERY_PATH)
print(f"  {query.n_obs} cells x {query.n_vars} genes")
print(f"  Layers available: {list(query.layers.keys())}")
if "counts" in query.layers:
    print(f"  'counts' layer max value: {query.layers['counts'].max()}")
else:
    print("  WARNING: no 'counts' layer found in query data either - "
          "we may need to check where raw counts are stored, or whether "
          "they need to be re-derived/re-exported from an earlier step.")

print("\n" + "=" * 70)
print("5. Shared genes between reference and query")
print("=" * 70)
shared_genes = set(ref.var_names) & set(query.var_names)
print(f"  Guilliams genes: {ref.n_vars}")
print(f"  Xiong query genes: {query.n_vars}")
print(f"  Shared genes: {len(shared_genes)}")
print(f"  ({100 * len(shared_genes) / query.n_vars:.1f}% of query genes are covered)")

print("\n" + "=" * 70)
print("SUMMARY - read before proceeding to training")
print("=" * 70)
print("If any WARNING appeared above, that needs to be resolved before "
      "writing the training script. Otherwise, paste this full output "
      "back and we'll design the training/inference scripts around your "
      "actual environment (CPU vs GPU affects epoch counts and expected "
      "runtime significantly).")