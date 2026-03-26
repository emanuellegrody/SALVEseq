# packages.R
# Installing and loading all required packages with version tracking


packages <- c(
  "clustree",
  "data.table",
  "doParallel",
  "dplyr",
  "foreach",
  "ggplot2",
  "ggrepel",
  "grid",
  "gridExtra",
  "metap",
  "msigdbr",
  "openxlsx",
  "parallel",
  "patchwork",
  "philentropy",
  "reshape2",
  "reticulate",
  "R.utils",
  "scales",
  "scCustomize",
  "Seurat",
  "SeuratWrappers",
  "stringdist",
  "stringr",
  "svglite",
  "tidyr",
  "transport"
)

versions_file <- "/gpfs/home/egy2296/SALVEseq/extractionScripts/R/packages_versions.csv"

load_packages_with_versions <- function(packages, versions_file, update_versions_file = TRUE) {
  
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", quiet = TRUE)
  }
  
  # Load existing versions file
  if (file.exists(versions_file)) {
    versions_df <- read.csv(versions_file, stringsAsFactors = FALSE)
  } else {
    versions_df <- data.frame(package = character(), version = character())
  }
  
  # ---- Phase 1: Install version-mismatched packages in a subprocess --------
  # A subprocess has a clean namespace environment, so unloadNamespace() will
  # never fail due to reverse-dependency imports from the main session.
  
  pkgs_to_install <- character()
  vers_to_install <- character()
  
  for (pkg in packages) {
    target_version <- versions_df$version[versions_df$package == pkg]
    if (length(target_version) == 0 || is.na(target_version)) next
    
    installed_version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) NA_character_
    )
    
    if (is.na(installed_version) || installed_version != target_version) {
      pkgs_to_install <- c(pkgs_to_install, pkg)
      vers_to_install <- c(vers_to_install, target_version)
    }
  }
  
  if (length(pkgs_to_install) > 0) {
    cat("Installing", length(pkgs_to_install),
        "package(s) with version mismatches...\n")
    
    # Each package is installed in its OWN Rscript process. This prevents
    # namespace contamination: installing CoGAPS in subprocess A may load
    # fgsea/data.table as transitive dependencies during dependency resolution.
    # If data.table were then installed in the same subprocess, the loaded
    # fgsea namespace would block unloadNamespace("data.table"). Separate
    # processes guarantee a clean namespace environment per install.
    #
    # dependencies = FALSE because these packages are already installed (just
    # at the wrong version). Their dependency trees are already satisfied on
    # disk. Using TRUE would trigger dependency resolution that can load
    # namespaces and recreate the original conflict.
    
    for (i in seq_along(pkgs_to_install)) {
      pkg <- pkgs_to_install[i]
      ver <- vers_to_install[i]
      cat("  [", i, "/", length(pkgs_to_install), "] ", pkg, " ", ver, " ...\n", sep = "")
      
      install_expr <- sprintf(
        'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", quiet = TRUE); remotes::install_version("%s", version = "%s", dependencies = FALSE, quiet = TRUE, upgrade = "never")',
        pkg, ver
      )
      
      exit_code <- system2("Rscript", args = c("-e", shQuote(install_expr)),
                           stdout = "", stderr = "")
      if (exit_code != 0) {
        cat("    WARNING: install of", pkg, ver, "exited with code", exit_code, "\n")
      }
    }
  }
  
  # ---- Phase 2: Load all packages and track new ones -----------------------
  
  new_packages <- character()
  new_versions <- character()
  
  for (pkg in packages) {
    target_version <- versions_df$version[versions_df$package == pkg]
    
    if (length(target_version) == 0 || is.na(target_version)) {
      # Package not in versions file: install latest if missing, then track
      installed_version <- tryCatch(
        as.character(packageVersion(pkg)),
        error = function(e) NA_character_
      )
      
      if (is.na(installed_version)) {
        cat("Installing", pkg, "(not in versions file, using latest)...\n")
        tryCatch({
          install.packages(pkg, quiet = TRUE, dependencies = TRUE)
          installed_version <- as.character(packageVersion(pkg))
        }, error = function(e) {
          cat("  Error installing", pkg, ":", e$message, "\n")
        })
      }
      
      if (!is.na(installed_version)) {
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        new_packages <- c(new_packages, pkg)
        new_versions <- c(new_versions, installed_version)
      }
      
    } else {
      # Package with pinned version: verify and load
      installed_version <- tryCatch(
        as.character(packageVersion(pkg)),
        error = function(e) NA_character_
      )
      
      if (!is.na(installed_version) && installed_version == target_version) {
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      } else {
        cat("WARNING:", pkg, "version", target_version, "not available after install.",
            "Installed:", ifelse(is.na(installed_version), "none", installed_version), "\n")
      }
    }
  }
  
  # ---- Update versions file with newly tracked packages --------------------
  
  if (update_versions_file && length(new_packages) > 0) {
    new_entries <- data.frame(package = new_packages, version = new_versions)
    versions_df <- rbind(versions_df, new_entries)
    versions_df <- versions_df[order(versions_df$package), ]
    write.csv(versions_df, versions_file, row.names = FALSE)
  }
  
  invisible(versions_df)
}

# Bioconductor packages (not available on CRAN)
bioc_packages <- c("fgsea", "zellkonverter", "SingleCellExperiment")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}

# Run the loader
load_packages_with_versions(packages, versions_file, update_versions_file = TRUE)

# Resolve namespace conflicts — ensure dplyr wins over MASS and plyr
if ("dplyr" %in% loadedNamespaces()) {
  conflict_fns <- c("select", "filter", "count", "rename", "mutate", "arrange", "summarise", "summarize")
  for (fn in conflict_fns) {
    if (exists(fn, envir = asNamespace("dplyr")))
      assign(fn, get(fn, envir = asNamespace("dplyr")), envir = .GlobalEnv)
  }
}

# ---- Python dependency: leidenalg (required for Seurat FindClusters algorithm = 4) ----
# Seurat calls leidenalg via reticulate. If the Python package is missing,
# FindClusters(..., algorithm = 4) will error at runtime.
#
# reticulate::py_install() installs into whichever Python reticulate has bound.
# On HPC systems, set RETICULATE_PYTHON in .Renviron or call
# reticulate::use_python("/path/to/python", required = TRUE) BEFORE this block
# to ensure the correct interpreter is targeted.

# leidenalg is optional — only needed for Seurat FindClusters(algorithm = 4).
# To install manually: pip install leidenalg