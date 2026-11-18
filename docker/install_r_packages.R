# install_r_packages.R
# Installs every R/Bioconductor package this project's runners depend on.
# Run once during the Docker build (see Dockerfile).

install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c(
    "SingleR",
    "celldex",
    "SingleCellExperiment",
    "scuttle",
    "BiocParallel",
    "BiocFileCache",
    "BiocSingular",
    "ExperimentHub",
    "AnnotationHub",
    "Matrix",
    "scrapper"   # needed by SingleR's aggregateReference() for large references
), update = FALSE, ask = FALSE)

install.packages(c("optparse", "dbplyr", "dplyr"),
                  repos = "https://cloud.r-project.org")

# SCINA was archived from CRAN (May 2026, administrative reason - see
# METHODS.md / project history) - install from the CRAN archive directly,
# NOT from GitHub HEAD, which has an unrelated bug (a hardcoded setwd()
# call not present in the actual released CRAN version).
#
# IMPORTANT: the exact version number below is NOT independently verified
# against what was actually used during development (we confirmed the
# install worked at the time, but never recorded the specific filename
# shown on the archive page). CHECK
# https://cran.r-project.org/src/contrib/Archive/SCINA/ and replace
# "SCINA_X.Y.Z.tar.gz" below with whatever version is actually listed
# there before building this image.
scina_version <- "SCINA_X.Y.Z.tar.gz"  # <- REPLACE with actual version from CRAN archive
scina_url <- paste0("https://cran.r-project.org/src/contrib/Archive/SCINA/", scina_version)
download.file(scina_url, destfile = scina_version)
install.packages(scina_version, repos = NULL, type = "source")

cat("\nAll R packages installed successfully.\n")
