"""
r_runner.py
Bridges Python (this wrapper) to R-based methods like SingleR by invoking
the R script as a subprocess and reading back its CSV output. This keeps
the wrapper genuinely multi-language rather than pretending SingleR can
run in pure Python.

RSCRIPT PATH: configurable via the RSCRIPT_BIN environment variable,
defaulting to this project's current Windows R install (R-4.6.1) if unset
- so existing Windows usage is unaffected. Inside the Docker image,
RSCRIPT_BIN is set to just "Rscript" (resolved via PATH, which the
container sets up correctly via apt) - see Dockerfile.
"""

import subprocess
import tempfile
import os
import pandas as pd
import scipy.io
import scipy.sparse as sp

DEFAULT_RSCRIPT_BIN = r"C:/Program Files/R/R-4.6.1/bin/Rscript.exe"


def run_r_method(adata, script_path: str, reference: str,
                  resolution: str = "main",
                  rscript_bin: str = None):
    """
    Exports adata's expression matrix + barcodes/genes to a temp directory,
    invokes the given R script as a subprocess, and reads back the
    resulting (barcode, label, confidence) CSV.

    Returns
    -------
    labels : pd.Series (indexed by adata.obs_names)
    confidence : pd.Series (indexed by adata.obs_names)
    """
    if rscript_bin is None:
        rscript_bin = os.environ.get("RSCRIPT_BIN", DEFAULT_RSCRIPT_BIN)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Export query data in the format run_singler.R expects
        scipy.io.mmwrite(os.path.join(tmp_dir, "query_counts.mtx"),
                          sp.csr_matrix(adata.X))
        pd.Series(adata.obs_names, name="barcode").to_csv(
            os.path.join(tmp_dir, "query_barcodes.csv"), index=False)
        pd.Series(adata.var_names, name="gene").to_csv(
            os.path.join(tmp_dir, "query_genes.csv"), index=False)

        output_csv = os.path.join(tmp_dir, "result.csv")

        cmd = [
            rscript_bin, script_path,
            "--reference", reference,
            "--resolution", resolution,
            "--input_dir", tmp_dir,
            "--output_csv", output_csv,
        ]
        print(f"  Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print(result.stderr)
            raise RuntimeError(f"R script failed (exit {result.returncode}): {script_path}")

        out_df = pd.read_csv(output_csv).set_index("barcode")

    labels = out_df["label"].reindex(adata.obs_names)
    confidence = out_df["conf"].reindex(adata.obs_names)
    return labels, confidence
