# Portability notes - read before building/running the Docker image

The Dockerfile sets up a reproducible Python + R environment, but several
scripts in this project have **hardcoded Windows-specific absolute
paths** from development on a specific machine. These will NOT work
inside the Linux container as-is and need to be changed first.

## Must fix before running in Docker

1. **`runners/scina_runner.py`** - `RSCRIPT_BIN` is hardcoded to
   `"C:/Program Files/R/R-4.6.1/bin/Rscript.exe"`. Change to:
   ```python
   RSCRIPT_BIN = "Rscript"
   ```
   (Inside the container, R is installed via apt and `Rscript` is
   already on PATH - no absolute path needed.)

2. **`runners/r_runner.py`** - check for a similarly hardcoded
   `RSCRIPT_BIN` path and apply the same fix.

3. **`config.py`** - several `MethodConfig` entries reference absolute
   Windows paths for trained model/reference files, e.g.:
   ```python
   kwargs={"model_name": r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\..."}
   ```
   These need to point at wherever the corresponding files are mounted
   inside the container instead (see "Data/model files" below).

## Data and model files

This repository's `.gitignore` deliberately excludes large files (`.h5ad`,
`.pt`, `scanvi_model/`, etc.) - they are NOT in version control and
therefore not in the Docker image either. To run the pipeline in a
container, mount a local directory containing these files as a volume:

```bash
docker run -v /path/to/your/data:/app/data your-image-name \
    python pipeline.py --input /app/data/your_file.h5ad --output /app/data/output.h5ad
```

The paths in `config.py` (model files, reference .h5ad files) should be
updated to point at this mounted `/app/data/...` location, or better,
made configurable via environment variables rather than hardcoded at
all - this is a reasonable next refactor once Docker use becomes routine
rather than occasional.

## Status

As of this writing, the Dockerfile builds a correct environment but the
path fixes above have NOT yet been applied to the actual runner/config
files. Treat this as infrastructure-in-progress, not a drop-in solution.
