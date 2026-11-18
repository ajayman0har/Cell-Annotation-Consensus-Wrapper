# Dockerfile - Cell Annotation Consensus Wrapper
#
# Builds a reproducible environment with both the Python stack
# (scanpy/scvi-tools/celltypist) and R/Bioconductor stack (SingleR/SCINA)
# this project depends on. Linux-based, so the painful Windows
# compilation issues encountered during development (annoy, tiledbsoma,
# harmonypy needing MSVC; SCINA needing a manual CRAN-archive install)
# mostly don't apply here - Linux has proper binary/source support for
# all of these via standard build tools.

FROM python:3.12-slim

# ── System dependencies + R ─────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libhdf5-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ── Python dependencies ──────────────────────────────────────────────────
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ── R / Bioconductor dependencies ───────────────────────────────────────
COPY docker/install_r_packages.R .
RUN Rscript install_r_packages.R

# ── Project code ─────────────────────────────────────────────────────────
COPY . .

# These override the Windows-default paths in r_runner.py/scina_runner.py
# (RSCRIPT_BIN) and config.py (ANNOTATION_WRAPPER_DATA_ROOT) - see those
# files' comments for the defaults they fall back to outside Docker.
ENV RSCRIPT_BIN=Rscript
ENV ANNOTATION_WRAPPER_DATA_ROOT=/app/data

CMD ["python", "pipeline.py", "--help"]
