# Portability notes - read before building/running the Docker image

The Dockerfile sets up a reproducible Python + R environment. The
previously hardcoded Windows-specific absolute paths have been made
environment-variable configurable (defaulting to the original Windows
paths, so existing Windows usage is unaffected), and the Dockerfile sets
the overrides needed for Linux:

1. **`runners/r_runner.py` and `runners/scina_runner.py`** - `RSCRIPT_BIN`
   now reads from the `RSCRIPT_BIN` environment variable, defaulting to
   the Windows R-4.6.1 path if unset. The Dockerfile sets
   `ENV RSCRIPT_BIN=Rscript`, which resolves correctly via PATH since R
   is installed via apt in the container.

2. **`config.py`** - model/reference file paths are now built from a
   `DATA_ROOT` variable, itself read from the
   `ANNOTATION_WRAPPER_DATA_ROOT` environment variable (defaulting to the
   Windows project path if unset). The Dockerfile sets this to
   `/app/data` - see "Data/model files" below for what needs to be
   mounted there.

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

Path portability fixes are applied (see above) and the Dockerfile sets
the needed environment variable overrides. Not yet verified by an actual
Docker build/run on this project's real data - the SCINA version number
in `docker/install_r_packages.R` still needs confirming against the live
CRAN archive listing (see that file's comment) before a build will
succeed end to end.
