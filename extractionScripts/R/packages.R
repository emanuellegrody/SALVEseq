# packages.R
# Installing and loading all required packages with version tracking

packages <- c(
  "CoGAPS",
  "data.table",
  "doParallel",
  "dplyr",
  "foreach",
  "ggplot2",
  "ggrepel",
  "grid",
  "gridExtra",
  "metap",
  "openxlsx",
  "parallel",
  "patchwork",
  "philentropy",
  "reshape2",
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

versions_file <- "/gpfs/home/egy2296/SALVEseq/packages_versions.csv"

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
  
  new_packages <- character()
  new_versions <- character()
  
  for (pkg in packages) {
    target_version <- versions_df$version[versions_df$package == pkg]
    
    if (length(target_version) == 0 || is.na(target_version)) {
      # Package not in versions file: install latest and track
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
      # Package in versions file: enforce specific version
      installed_version <- tryCatch(
        as.character(packageVersion(pkg)),
        error = function(e) NA_character_
      )
      
      if (!is.na(installed_version) && installed_version == target_version) {
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      } else {
        tryCatch({
          remotes::install_version(
            pkg,
            version = target_version,
            dependencies = TRUE,
            quiet = TRUE,
            upgrade = "never"
          )
          library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        }, error = function(e) {
          cat("  Error installing", pkg, "version", target_version, ":", e$message, "\n")
        })
      }
    }
  }
  
  # Update versions file with new packages
  if (update_versions_file && length(new_packages) > 0) {
    new_entries <- data.frame(package = new_packages, version = new_versions)
    versions_df <- rbind(versions_df, new_entries)
    versions_df <- versions_df[order(versions_df$package), ]
    write.csv(versions_df, versions_file, row.names = FALSE)
  }
  
  invisible(versions_df)
}

# Run the loader
load_packages_with_versions(packages, versions_file, update_versions_file = TRUE)