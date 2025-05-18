# setup.R
# installing and loading all required packages

# CRAN
packages <- c(
  "stringr",
  "data.table",
  "dplyr",
  "tidyr",
  "patchwork",
  "ggplot2",
  "scales",
  "svglite",
  "Seurat",
  "scCustomize",
  "metap",
  "reshape2",
  "openxlsx",
  "R.utils",
  "stringdist",
  "parallel",
  "doParallel",
  "foreach",
  "openxlsx"
)


for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

