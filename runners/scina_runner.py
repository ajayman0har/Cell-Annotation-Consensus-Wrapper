"""
scina_runner.py
Bridges Python (this wrapper) to the SCINA R script, mirroring
r_runner.py's export/subprocess pattern but with a simpler CLI - SCINA
needs no "reference"/"resolution" concept, just the query data and its
own baked-in marker signature panel (see run_scina.R).
"""

import subprocess
import tempfile
import os
import pandas as pd
import scipy.io
import scipy.sparse as sp

# Same R install used for SingleR - see r_runner.py's note for why this
# is hardcoded rather than relying on PATH.
RSCRIPT_BIN = r"C:/Program Files/R/R-4.6.1/bin/Rscript.exe"


def run_scina(adata, script_path: str = "runners/run_scina.R"):
    """
    Returns
    -------
    labels : pd.Series (indexed by adata.obs_names)
    confidence : pd.Series (indexed by adata.obs_names)
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        scipy.io.mmwrite(os.path.join(tmp_dir, "query_counts.mtx"),
                          sp.csr_matrix(adata.X))
        pd.Series(adata.obs_names, name="barcode").to_csv(
            os.path.join(tmp_dir, "query_barcodes.csv"), index=False)
        pd.Series(adata.var_names, name="gene").to_csv(
            os.path.join(tmp_dir, "query_genes.csv"), index=False)

        output_csv = os.path.join(tmp_dir, "result.csv")

        cmd = [
            RSCRIPT_BIN, script_path,
            "--input_dir", tmp_dir,
            "--output_csv", output_csv,
        ]
        print(f"  Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print(result.stderr)
            raise RuntimeError(f"SCINA R script failed (exit {result.returncode})")

        out_df = pd.read_csv(output_csv).set_index("barcode")

    labels = out_df["label"].reindex(adata.obs_names)
    confidence = out_df["conf"].reindex(adata.obs_names)
    return labels, confidence