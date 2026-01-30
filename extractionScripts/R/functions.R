# functions.R

# Loading and Processing
## SALVE and GEX: extract category from filename
extract_categoryv6 <- function(filename) {
  # Split the filename into parts
  target_part <- basename(filename)
  target_part <- strsplit(target_part, "_bamsort_alignment_")[[1]][1]
  alignment_part <- strsplit(filename, "_bamsort_alignment_")[[1]][2]
  alignment_part <- sub("\\.csv$", "", alignment_part)
  
  targets <- strsplit(target_part, "_")[[1]]
  targets <- targets[-1]
  targets <- targets[!targets %in% c("PL", "CI")]
  if (all(is.na(targets))) {
    return(NA)
  }
  
  # universal truths
  if (grepl("KLRB1$", alignment_part)) return("KLRB1")
  if (length(targets) == 1 & "pol" %in% targets) return("US")
  if (grepl("LTR_D1$", alignment_part)) return("any-D1")
  if (grepl("D1_A1$", alignment_part)) return("US")
  
  
  # treat GEX like SALVE env/SSenv/pol
  if (alignment_part == "A1_D4") {
    if ("D1" %in% targets) return("SS-MS4")
    else return("any-D4")
  }
  if (alignment_part == "D4_A7") {
    if ("D1" %in% targets) return("SS")
    else return("US-SS")
  }
  if (alignment_part == "A7_LTR") {
    if ("D1" %in% targets) return("SS-MS7")
    if ("SSenv" %in% targets) return("US-SS")
    if ("tat" %in% targets | "D4" %in% targets) return("MS")
    else return("any-A7")
  }
  
  # Default case
  cat("  Failed to categorize. Assigning as: NA\n")
  return(NA)
}
## SALVE: spliced reads join 
## Load correction CSV and apply category overrides
apply_umi_corrections <- function(umi_read_counts, correction_file, target_category = "MS") {
  # Check if correction file exists
  if (!file.exists(correction_file)) {
    cat("No correction file found at:", correction_file, "\n")
    return(umi_read_counts)
  }
  
  allowed_categories <- c("SS-MS4", "SS-MS7", "any-A7", "any-D4")
  # Load correction data
  corrections <- read.csv(correction_file)
  
  # Validate required columns
  required_cols <- c("CB", "UB")
  if (!all(required_cols %in% colnames(corrections))) {
    cat("Warning: Correction file missing required columns (CB, UB)\n")
    return(umi_read_counts)
  }
  
  # Create composite key for matching
  corrections$key <- paste(corrections$CB, corrections$UB, sep = "_")
  umi_read_counts$key <- paste(umi_read_counts$cellID, umi_read_counts$UMI, sep = "_")
  
  in_correction_file <- umi_read_counts$key %in% corrections$key
  
  # Find matches with allowed categories
  to_correct <- in_correction_file & 
    umi_read_counts$category %in% allowed_categories
  to_toss <- in_correction_file & 
    !umi_read_counts$category %in% allowed_categories
  
  n_corrected <- sum(to_correct)
  n_tossed <- sum(to_toss)
  
  
  if (n_corrected > 0) {
    # Store original categories for reporting
    original_categories <- umi_read_counts$category[to_correct]
    
    cat("Corrected to", target_category, n_corrected, "UMI from:", 
        paste(unique(original_categories), collapse = ", "), "\n")
  } else {
    cat("No matching UMIs found in correction file\n")
  }
  
  if (n_tossed > 0) {
    tossed_categories <- umi_read_counts$category[to_toss]
    cat("Tossed", n_tossed, "UMIs with non-allowed categories:", 
        paste(unique(tossed_categories), collapse = ", "), "\n")
  }
  
  # Apply correction
  umi_read_counts$category[to_correct] <- target_category
  umi_read_counts <- umi_read_counts[!to_toss, ]
  
  # Remove temporary key column
  umi_read_counts$key <- NULL
  
  return(umi_read_counts)
}
## SALVE and GEX: resolve UMI mapping to multiple genes/regions
resolve_multimapv6 <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  valid_pairs <- list(
    # keep the second entry
    c("any-D1", "US"),
    c("any-D4", "US"),
    c("SS-MS4", "US"),
    c("any-D4", "US-SS"), 
    c("any-A7", "US-SS"), 
    c("MS", "US-SS"), 
    c("SS-MS4", "SS"),
    c("SS-MS7", "SS"),
    c("any-D4", "SS-MS4"),
    c("any-A7", "SS-MS7"),
    c("US-SS", "SS"),
    c("SS-MS7", "MS"),
    c("any-A7", "MS")
  )
  
  # Create a composite key
  df$key <- paste(df$cellID, df$UMI, sep = "_")
  
  # Identify multi-category keys
  key_counts <- table(df$key)
  multi_keys <- names(key_counts[key_counts > 1])
  
  # Statistics counters
  total_multimapped <- length(multi_keys)
  resolved_count <- 0
  tossed_count <- 0
  tossed_more_than_2_count <- 0
  consolidated_same_category <- 0 
  
  # Split into single and multi entries
  singles <- df[!df$key %in% multi_keys, ]
  multis <- df[df$key %in% multi_keys, ]
  
  # Initialize result as empty data frame with correct structure
  result <- data.frame(
    cellID = character(0),
    UMI = character(0),
    category = character(0),
    read_count = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Add singles to result
  if (nrow(singles) > 0) {
    singles_clean <- data.frame(
      cellID = as.character(singles$cellID),
      UMI = as.character(singles$UMI),
      category = as.character(singles$category),
      read_count = as.numeric(singles$read_count),
      stringsAsFactors = FALSE
    )
    result <- rbind(result, singles_clean)
  }
  
  # Process each unique key in multis
  if (nrow(multis) > 0) {
    unique_keys <- unique(multis$key)
    
    for (k in unique_keys) {
      # Get all rows for this key
      rows <- multis[multis$key == k, ]
      categories <- unique(rows$category)
      
      # Process based on number of categories
      if (length(categories) == 1) {
        # For keys with 1 category but multiple rows, consolidate by summing read_counts
        total_reads <- sum(as.numeric(rows$read_count))
        keep_row <- data.frame(
          cellID = as.character(rows$cellID[1]),
          UMI = as.character(rows$UMI[1]),
          category = as.character(rows$category[1]),
          read_count = total_reads,
          stringsAsFactors = FALSE
        )
        result <- rbind(result, keep_row)
        consolidated_same_category <- consolidated_same_category + 1
      } else if (length(categories) == 2) {
        # Check against valid pairs
        valid_match <- FALSE
        for (pair in valid_pairs) {
          if (all(sort(categories) == sort(pair))) {
            # Valid pair, keep the second category
            keep_cat <- pair[2]
            keep_rows <- rows[rows$category == keep_cat, ]
            keep_row <- data.frame(
              cellID = as.character(keep_rows$cellID[1]),
              UMI = as.character(keep_rows$UMI[1]),
              category = as.character(keep_cat),
              read_count = sum(as.numeric(rows$read_count)),
              stringsAsFactors = FALSE
            )
            result <- rbind(result, keep_row)
            resolved_count <- resolved_count + 1
            valid_match <- TRUE
            break
          }
        }
        
        # If no valid match found, count as tossed
        if (!valid_match) {
          tossed_count <- tossed_count + 1
        }
      } else if (length(categories) > 2) {
        # Toss UMIs with more than 2 categories
        tossed_more_than_2_count <- tossed_more_than_2_count + 1
      }
    }
  }
  
  # Print statistics
  cat("Total multi-mapped UMIs:", total_multimapped, "\n")
  cat("UMIs successfully resolved:", resolved_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, resolved_count/total_multimapped*100, 0)), "\n")
  cat("UMIs tossed (no valid pair):", tossed_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, tossed_count/total_multimapped*100, 0)), "\n")
  cat("UMIs tossed (>2 categories):", tossed_more_than_2_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, tossed_more_than_2_count/total_multimapped*100, 0)), "\n")
  cat("Total UMIs tossed:", tossed_count + tossed_more_than_2_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, 
                                 (tossed_count + tossed_more_than_2_count)/total_multimapped*100, 0)), "\n")
  
  # Clean up row names
  rownames(result) <- NULL
  
  return(result)
}
## SALVE and GEX: read and process BAM
process_bamsortv6 <- function(samples_list, input.dir, output.dir, raw_cellIDs, 
                              correction_dir = NULL) {
  #
  # optional to read in paste0(sample, "_D1_splicesites_filtered.csv")
  # which is output from bamsort_splice_sites to update MS
  #
  if (class(samples_list) != "character") {
    stop("samples_list input must be a list of sample names")
  }
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  for (sample in samples_list) {
    cat("\nProcessing sample:", sample, "\n")
    
    # Get all files for current sample
    sample_files <- unlist(lapply(input.dir, function(dir) {
      list.files(dir, 
                 pattern = paste0(".*", sample, ".*_bamsort_alignment_.*\\.csv$"),
                 full.names = TRUE)
    }))
    
    if (length(sample_files) == 0) {
      cat("No files found for sample:", sample, "\n")
      next
    }
    
    # Extract categories for each file
    categories <- sapply(sample_files, extract_categoryv6)
    sample_files <- sample_files[!is.na(categories)]
    
    # Create a data frame to store the combined results
    all_data <- data.frame()
    filenames <- basename(sample_files)
    
    # Process each relevant file
    for (j in seq_along(sample_files)) {
      extracted_data <- NULL
      
      tryCatch({
        file_data <- fread(sample_files[j], data.table = FALSE)
        
        if(nrow(file_data) == 0) {
          cat("Warning: File is empty:", filenames[j], "\n")
          next
        }
        
        required_cols <- c("cellID", "UMI", "count")
        missing_cols <- setdiff(required_cols, colnames(file_data))
        
        if (length(missing_cols) > 0) {
          cat("Warning: Missing required columns:", paste(missing_cols, collapse=", "), "\n")
          cat("Available columns:", paste(colnames(file_data), collapse=", "), "\n")
          next
        }
        
        extracted_data <- file_data %>%
          select(cellID, UMI, count) %>%
          rename(read = count) %>%
          mutate(category = categories[j])
        
      }, error = function(e) {
        cat("Error reading file:", filenames[j], "\nLikely bad data file\n")
        cat("Error message:", conditionMessage(e), "\n")
      })
      
      if (!is.null(extracted_data) && nrow(extracted_data) > 0) {
        all_data <- rbind(all_data, extracted_data)
      }
    }
    
    # If we have data, process it
    if (nrow(all_data) > 0) {
      all_data <- all_data %>% 
        filter(!is.na(cellID) & !is.na(UMI) & !is.na(read)) %>%
        unique()
      
      tryCatch({
        umi_read_counts <- all_data %>%
          group_by(cellID, UMI, category) %>%
          summarize(
            read_count = sum(as.numeric(read)),
            .groups = 'drop'
          )
        
        cat("Gathered read counts for", nrow(umi_read_counts), "UMIs from", 
            n_distinct(umi_read_counts$cellID), "cells\n")
        
        # Create a wide format with UMI counts per category
        umi_by_category <- umi_read_counts %>%
          group_by(cellID, category) %>%
          summarize(
            category_UMIs = n_distinct(UMI),
            category_reads = sum(read_count),
            .groups = 'drop'
          )
        
        umi_by_category_wide <- tidyr::pivot_wider(
          umi_by_category,
          id_cols = cellID,
          names_from = category,
          names_sep = "_",
          values_from = c(category_UMIs, category_reads),
          values_fill = 0)
        
        # Apply raw_cellIDs filtering
        if (sample %in% names(raw_cellIDs) && !is.null(raw_cellIDs[[sample]])) {
          valid_cells <- raw_cellIDs[[sample]]
          cat("Keeping only cells in", 
              length(valid_cells), "valid cells from 10X data\n")
          
          valid_umi_read_counts <- umi_read_counts %>%
            filter(cellID %in% valid_cells)
          
          # Clean category values
          if (!is.null(names(valid_umi_read_counts$category))) {
            category_values <- as.character(valid_umi_read_counts$category)
            category_values <- gsub('^"(.*)"$', '\\1', category_values)
            valid_umi_read_counts$category <- category_values
          } else {
            valid_umi_read_counts$category <- gsub('^"(.*)"$', '\\1', as.character(valid_umi_read_counts$category))
          }
          
          # Look for splice correction file
          if (!is.null(correction_dir)) {
            correction_file <- file.path(correction_dir, paste0(sample, "_D1_splicesites_filtered.csv"))
            valid_umi_read_counts <- apply_umi_corrections(
              valid_umi_read_counts, 
              correction_file, 
              target_category = "MS"
            )
          }
          
          # Multimap resolution (with corrected categories)
          filtered_umi_read_counts <- resolve_multimapv6(valid_umi_read_counts)
          
          # Convert to basic types
          filtered_umi_read_counts <- data.frame(
            cellID = as.character(filtered_umi_read_counts$cellID),
            UMI = as.character(filtered_umi_read_counts$UMI),
            category = as.character(filtered_umi_read_counts$category),
            read_count = as.numeric(filtered_umi_read_counts$read_count),
            stringsAsFactors = FALSE
          )
          
          # Consolidate any categories
          filtered_umi_read_counts <- filtered_umi_read_counts %>%
            mutate(category = if_else(grepl("^any", category), "any", category)) %>%
            mutate(category = if_else(grepl("^SS-MS", category), "SS-MS", category))
          
          cat("After filtering: kept", nrow(filtered_umi_read_counts), "UMIs from", 
              n_distinct(filtered_umi_read_counts$cellID), "valid cells\n")
        } else {
          cat("Warning: No 10X data found for sample", sample, "- using unfiltered cell list\n")
          filtered_umi_read_counts <- umi_read_counts
        }
        
        # Convert to basic types
        all_data_clean <- data.frame(
          cellID = as.character(all_data$cellID),
          UMI = as.character(all_data$UMI),
          read = as.numeric(all_data$read),
          category = as.character(all_data$category),
          stringsAsFactors = FALSE
        )
        
        umi_read_counts_clean <- data.frame(
          cellID = as.character(umi_read_counts$cellID),
          UMI = as.character(umi_read_counts$UMI),
          category = as.character(umi_read_counts$category),
          read_count = as.numeric(umi_read_counts$read_count),
          stringsAsFactors = FALSE
        )
        
        # Save outputs
        umi_level_file <- file.path(output.dir, paste0(sample, "_UMI_read_counts_raw.csv"))
        write.csv(umi_read_counts_clean, file = umi_level_file, row.names = FALSE)
        
        filtered_umi_level_file <- file.path(output.dir, paste0(sample, "_UMI_read_counts_full.csv"))
        write.csv(filtered_umi_read_counts, file = filtered_umi_level_file, row.names = FALSE)
        
        cat("Successfully processed\n")
        rm(umi_by_category, umi_by_category_wide, umi_read_counts)
      }, error = function(e) {
        cat("Error processing data for sample", sample, ":", conditionMessage(e), "\n")
        
        tryCatch({
          all_data_clean <- data.frame(
            cellID = as.character(all_data$cellID),
            UMI = as.character(all_data$UMI),
            read = as.numeric(all_data$read),
            category = as.character(all_data$category),
            stringsAsFactors = FALSE
          )
          
          raw_file <- file.path(output.dir, paste0(sample, "_raw_data.csv"))
          write.csv(all_data_clean, file = raw_file, row.names = FALSE)
          cat("Saved raw data to:", raw_file, "\n")
        }, error = function(e2) {
          cat("Could not save raw data:", conditionMessage(e2), "\n")
        })
      })
    } else {
      cat("No data extracted for sample:", sample, "\n")
    }
  }
}
## SALVE and GEX: filtering data for minimums
set_minimumsv6 <- function(mode, umi_file, min_reads = 1, min_region_count = 1, 
                           min_umi = 1, min_cells = 0, min_reads_cell = 1,
                           save_valid_umi = TRUE, valid_umi_output_file = NULL) {
  
  # Define categories based on mode
  categories <- switch(mode,
                       "SALVE" = c("US", "US-SS", "SS-MS", "SS", "MS", "any"),
                       "GEX" = c("US", "US-SS", "any"),
                       "KLRB1" = c("KLRB1"),
                       stop("Invalid mode. Must be 'SALVE', 'GEX', or 'KLRB1'")
  )
  
  # Helper function to return empty result
  return_empty <- function(message = "", umi_entries = 0) {
    if (nzchar(message)) cat(message, "\n")
    empty_data <- data.frame(cellID = character(0), stringsAsFactors = FALSE)
    empty_data[categories] <- lapply(categories, function(x) numeric(0))
    empty_data$total <- numeric(0)
    cat("Converted", umi_entries, "UMI entries to 0 cells\n")
    return(empty_data)
  }
  
  # Read and validate data
  umi_data <- read.csv(umi_file, stringsAsFactors = FALSE)
  if (nrow(umi_data) == 0) {
    return(return_empty("No UMI data found in file"))
  }
  
  original_entries <- nrow(umi_data)
  
  # Apply read count filter and check for empty result
  if (min_reads > 0) {
    umi_data <- umi_data[umi_data$read_count >= min_reads, ]
    if (nrow(umi_data) == 0) {
      return(return_empty("All UMI data filtered out by min_reads threshold", original_entries))
    }
  }
  
  # Calculate cell-level metrics using dplyr for efficiency
  cell_stats <- umi_data %>%
    group_by(cellID) %>%
    summarize(
      total_reads = sum(read_count, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Apply reads per cell filter early
  if (min_reads_cell > 0) {
    valid_cells <- cell_stats$cellID[cell_stats$total_reads >= min_reads_cell]
    umi_data <- umi_data[umi_data$cellID %in% valid_cells, ]
    if (nrow(umi_data) == 0) {
      return(return_empty("All cells filtered out by min_reads_cell threshold", original_entries))
    }
  }
  
  # Create UMI count matrix more efficiently
  umi_counts <- umi_data %>%
    count(cellID, category, name = "umi_count") %>%
    pivot_wider(names_from = category, values_from = umi_count, values_fill = 0)
  
  # Initialize result with all categories (including missing ones as 0)
  result <- umi_counts %>%
    select(cellID, any_of(categories)) %>%
    mutate(across(-cellID, ~replace_na(.x, 0)))
  
  # Add missing categories as zero columns
  missing_cats <- setdiff(categories, colnames(result))
  if (length(missing_cats) > 0) {
    result[missing_cats] <- 0
  }
  
  # Reorder columns to match expected order
  result <- result[c("cellID", categories)]
  
  # Apply minimum cells per region filter
  if (min_cells > 0) {
    for (cat in categories) {
      if (sum(result[[cat]] > 0) < min_cells) {
        result[[cat]] <- 0
      }
    }
  }
  
  # Calculate metrics for filtering
  result <- result %>%
    mutate(
      region_count = rowSums(across(all_of(categories), ~ .x > 0)),
      total = rowSums(across(all_of(categories)))
    )
  
  # Apply filters
  result <- result %>%
    filter(
      region_count >= min_region_count,
      total >= min_umi
    ) %>%
    select(-region_count)  # Remove helper column
  
  # Check for empty result
  if (nrow(result) == 0) {
    return(return_empty("All cells filtered out by region_count or UMI thresholds", original_entries))
  }
  
  # Save valid cellID-UMI pairs if requested
  if (save_valid_umi && !is.null(valid_umi_output_file) && nrow(result) > 0) {
    # Extract valid cellIDs after all filtering
    valid_cellIDs <- result$cellID
    
    # Filter original UMI data to keep only valid cell-UMI combinations
    valid_umi_data <- umi_data %>%
      filter(cellID %in% valid_cellIDs) %>%
      select(cellID, UMI)
    
    # Write to CSV without row names for cleaner output
    write.csv(valid_umi_data, file = valid_umi_output_file, row.names = FALSE)
    cat("Saved", nrow(valid_umi_data), "valid cellID-UMI pairs\n")
  }
  
  cat("Converted", original_entries, "UMI entries to", nrow(result), "cells\n")
  return(result)
}
## SALVE and GEX: filter multiple samples
process_all_set_minimumsv6 <- function(mode, sample_list, input_dir, save_output = TRUE, 
                                       output_dir = "/projects/b1042/GoyalLab/egrody/extractedData/",
                                       min_reads = 1, min_region_count = 1, min_umi = 1, 
                                       min_cells = 0, min_reads_cell = 1,
                                       save_valid_umi = TRUE, valid_umi_dir = NULL,
                                       aggregate_valid_umi = TRUE) {
  
  # Create output directory if needed
  if (save_output == TRUE) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create separate directory for valid UMI files if needed
  # Use the provided output_dir as the base for valid_umi_dir
  if (save_valid_umi == TRUE) {
    if (is.null(valid_umi_dir)) {
      valid_umi_dir <- file.path(output_dir, "valid_UMI")
    }
    dir.create(valid_umi_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Write filters to file
  filters_file <- file.path(output_dir, "filters.txt")
  filter_params <- c(
    paste("min_reads =", min_reads),
    paste("min_region_count =", min_region_count),
    paste("min_umi =", min_umi),
    paste("min_cells =", min_cells),
    paste("min_reads_cell =", min_reads_cell),
    paste("mode =", mode),
    paste("processed_date =", Sys.Date())
  )
  writeLines(filter_params, filters_file)
  
  # Define output suffix based on mode
  suffix <- paste0("_", mode, "_filtered.csv")
  umi_suffix <- paste0("_", mode, "_valid_cellID_UMI.csv")
  
  # Initialize list to collect valid UMI data if aggregating
  all_valid_umi <- if (aggregate_valid_umi && save_valid_umi) list() else NULL
  
  # Process samples with better error handling
  results <- sapply(sample_list, function(sample) {
    cat("\nProcessing sample:", sample, "\n")
    
    input_file <- file.path(input_dir, paste0(sample, "_UMI_read_counts_full.csv"))
    output_file <- file.path(output_dir, paste0(sample, suffix))
    
    # Construct valid UMI output filename with sample name and mode included
    valid_umi_file <- if (save_valid_umi) {
      file.path(valid_umi_dir, paste0(sample, umi_suffix))
    } else {
      NULL
    }
    
    # Check if input exists
    if (!file.exists(input_file)) {
      warning("Input file not found: ", input_file)
      return("failed")
    }
    
    # Process with error handling
    tryCatch({
      result <- set_minimumsv6(
        mode = mode, 
        umi_file = input_file, 
        min_reads = min_reads, 
        min_region_count = min_region_count, 
        min_umi = min_umi, 
        min_cells = min_cells,
        min_reads_cell = min_reads_cell,
        save_valid_umi = save_valid_umi,
        valid_umi_output_file = valid_umi_file
      )
      
      # Save filtered cell counts
      if (save_output == TRUE) {
        write.csv(result, file = output_file, row.names = FALSE)
      }
      
      # Collect valid UMI data for aggregation
      if (aggregate_valid_umi && save_valid_umi && file.exists(valid_umi_file)) {
        sample_umi_data <- read.csv(valid_umi_file, stringsAsFactors = FALSE)
        sample_umi_data$sample <- sample
        all_valid_umi[[sample]] <<- sample_umi_data
      }
      
      return("success")
      
    }, error = function(e) {
      cat("Error processing sample", sample, ":", conditionMessage(e), "\n")
      return("failed")
    })
  }, USE.NAMES = FALSE)
  
  # Write aggregated valid UMI file if requested
  if (aggregate_valid_umi && save_valid_umi && length(all_valid_umi) > 0) {
    aggregated_file <- file.path(output_dir, paste0("all_samples_", mode, "_valid_cellID_UMI.csv"))
    aggregated_data <- do.call(rbind, all_valid_umi)
    
    # Reorder columns: sample, cellID, UMI
    aggregated_data <- aggregated_data[, c("sample", "cellID", "UMI")]
    
    write.csv(aggregated_data, file = aggregated_file, row.names = FALSE)
    cat("\nSaved aggregated valid UMI data with", nrow(aggregated_data), 
        "entries across", length(all_valid_umi), "samples\n")
  }
  
  # Summary
  successful <- sum(results == "success")
  failed <- sum(results == "failed")
  cat("\nProcessing complete! Successful:", successful, "Failed:", failed, "\n")
  
  if (save_output == TRUE) {
    return(invisible(list(successful = successful, failed = failed)))
  }
  else {
    return(result)
  }
}
## Splicing
classify_site <- function(value, reference_values, reference_sites) {
  if (is.na(value)) {
    return("nc")
  }
  
  distances <- abs(value - reference_values)
  valid_distances <- distances[!is.na(distances)]
  
  if (length(valid_distances) == 0) {
    return("nc")
  }
  
  min_dist <- min(valid_distances)
  closest_site <- reference_sites[which.min(distances)]
  
  if (min_dist <= 3) {
    return(as.character(closest_site))
  } else if (min_dist <= 12) {
    return("mnc")
  } else {
    return("nc")
  }
}

## Combined GEX and SALVE
resolve_multimap_gs <- function(df) {
  # Priority categories that override others
  priority_cats <- c("US", "SS", "MS")
  
  # Initialize counters
  stats <- list(
    total = 0,
    resolved_same = 0,
    resolved_priority = 0,
    resolved_to_any = 0,
    no_conflict = 0
  )
  
  # Process each row
  result <- vector("list", nrow(df))
  
  for (i in seq_len(nrow(df))) {
    cat_x <- df$category.x[i]
    cat_y <- df$category.y[i]
    # read_x <- ifelse(is.na(df$read_count.x[i]), 0, df$read_count.x[i])
    # read_y <- ifelse(is.na(df$read_count.y[i]), 0, df$read_count.y[i])
    
    # Case 1: Both categories are NA - skip this row
    if (is.na(cat_x) && is.na(cat_y)) {
      result[[i]] <- NULL
      next
    }
    
    # Case 2: One category is NA - use the non-NA category
    if (is.na(cat_x) && !is.na(cat_y)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_y,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$no_conflict <- stats$no_conflict + 1
      next
    }
    
    if (!is.na(cat_x) && is.na(cat_y)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_x,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$no_conflict <- stats$no_conflict + 1
      next
    }
    
    # Case 3: Both categories present - apply resolution rules
    stats$total <- stats$total + 1
    
    # Rule 0: If both categories are identical, keep that category
    if (cat_x == cat_y) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_x,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$resolved_same <- stats$resolved_same + 1
      next
    }
    
    # Rule 1: If one is a priority category (US, SS, MS), keep that one
    if (cat_x %in% priority_cats && !(cat_y %in% priority_cats)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_x,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$resolved_priority <- stats$resolved_priority + 1
      next
    }
    
    if (cat_y %in% priority_cats && !(cat_x %in% priority_cats)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_y,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$resolved_priority <- stats$resolved_priority + 1
      next
    }
    
    # Rule 2: Otherwise (both are priority categories or neither are), resolve to "any"
    result[[i]] <- data.frame(
      cellID = df$cellID[i],
      UMI = df$UMI[i],
      category = "any",
      #read_count = read_x + read_y,
      stringsAsFactors = FALSE
    )
    stats$resolved_to_any <- stats$resolved_to_any + 1
  }
  
  # Combine results, excluding NULLs
  result <- do.call(rbind, result[!sapply(result, is.null)])
  rownames(result) <- NULL
  
  # Print statistics
  cat("Multi-mapped UMIs (both categories present):", stats$total, "\n")
  cat("UMIs with identical categories:", stats$resolved_same, 
      sprintf("(%.1f%%)", ifelse(stats$total > 0, stats$resolved_same/stats$total*100, 0)), "\n")
  cat("UMIs resolved by priority category (US/SS/MS):", stats$resolved_priority, 
      sprintf("(%.1f%%)", ifelse(stats$total > 0, stats$resolved_priority/stats$total*100, 0)), "\n")
  cat("UMIs resolved to 'any':", stats$resolved_to_any, 
      sprintf("(%.1f%%)", ifelse(stats$total > 0, stats$resolved_to_any/stats$total*100, 0)), "\n")
  cat("UMIs with no conflict (one or both categories NA):", stats$no_conflict, "\n")
  
  
  return(result)
}


# Loading
## standard processing of counts matrices with Seurat
SeuratPipeline <- function(file_location, sample_name, output_dir = "/projects/b1042/GoyalLab/egrody/", 
                           plots = FALSE, rds = FALSE, species = "macaque", remove_mac = TRUE) {
  #if you would like to output plots from this pipeline, pass entries for the plots and output_dir variables
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Reading in feature matrices
  data <- Read10X(data.dir = file_location)
  # Initialize the Seurat object with the raw (non-normalized data)
  sample <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  # QC
  sample <- Add_Mito_Ribo(sample, species = species)
  if (plots) {
    vplot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
    ggsave(vplot, file = paste0(output_dir, sample_name, "_preQC_violinplot.svg"))
  }
  sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA < 20000 & percent_mito < 5)
  if (plots) {
    vplot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
    ggsave(vplot, file = paste0(output_dir, sample_name, "_postQC_violinplot.svg"))
  }
  # Remove CD4/CD8 doublets
  if ("CD4" %in% rownames(sample) && ("CD8A" %in% rownames(sample) || "CD8B" %in% rownames(sample))) {
    cd4_expr <- GetAssayData(sample, layer = "counts")["CD4", ]
    cd8a_expr <- if ("CD8A" %in% rownames(sample)) GetAssayData(sample, layer = "counts")["CD8A", ] else rep(0, ncol(sample))
    cd8b_expr <- if ("CD8B" %in% rownames(sample)) GetAssayData(sample, layer = "counts")["CD8B", ] else rep(0, ncol(sample))
    
    cd4_pos <- cd4_expr > 0
    cd8_pos <- (cd8a_expr > 0) | (cd8b_expr > 0)
    doublets <- cd4_pos & cd8_pos
    
    sample <- subset(sample, cells = colnames(sample)[!doublets])
  }
  # Remove mac239 gene
  if (remove_mac == TRUE) {
    if ("mac239" %in% rownames(sample)) {
      sample <- sample[!rownames(sample) %in% "mac239", ]
    }
  }
  # Normalizing
  sample <- NormalizeData(sample, verbose = FALSE)
  # Feature selection
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # Scaling
  genes <- rownames(sample)
  sample <- ScaleData(sample, features = genes, verbose = FALSE)
  # Linear dimension reduction
  sample <- RunPCA(sample, features = VariableFeatures(object = sample), verbose = FALSE)
  if (plots) {
    eplot <- ElbowPlot(sample, ndims = 50)
    ggsave(eplot, file = paste0(output_dir, sample_name, "_elbow.svg"))
  }
  # Clustering
  sample <- FindNeighbors(sample, dims = 1:30, verbose = FALSE)
  sample <- FindClusters(sample, resolution = 0.3, verbose = FALSE)
  # Dimension reduction
  sample <- RunUMAP(sample, dims = 1:30, verbose = FALSE)#, umap.method = "umap-learn", metric = "correlation")
  if (plots) {
    dplot <- DimPlot(sample, reduction = "umap")
    ggsave(dplot, file = paste0(output_dir, sample_name, "_umap.svg"))
  }
  if (rds) {
    saveRDS(sample, file = paste0(output_dir, sample_name, ".rds"))
  }
  
  # Return the variable
  return(sample)
}


# Output Dataframes
## isolate specific gene(s) expression from a Seurat object
targetExpressionDF <- function(seurat_obj, genes, count_type = "normalized") {
  # Ensure genes is a vector
  if (!is.vector(genes)) {
    genes <- c(genes)
  }
  
  # Get cell IDs
  cell_ids <- colnames(seurat_obj)
  
  # Initialize dataframe with cellID
  expression_df <- data.frame(
    cellID = cell_ids,
    stringsAsFactors = FALSE
  )
  
  # Check which genes are available
  available_genes <- rownames(seurat_obj)
  
  # Get expression data (choose raw or normalized)
  if (count_type == "raw") {
    expression_data <- GetAssayData(object = seurat_obj, assay = "RNA", layer = "counts")
  } else {
    expression_data <- GetAssayData(object = seurat_obj, assay = "RNA", layer = "data")
  }
  
  # Add column for each requested gene
  for (gene in genes) {
    if (gene %in% available_genes) {
      # Gene present: extract expression values
      expression_df[[gene]] <- as.numeric(expression_data[gene, ])
    } else {
      # Gene absent: add zeros and warn
      expression_df[[gene]] <- rep(0, length(cell_ids))
      warning(paste("Gene", gene, "not found in dataset. \n"))
    }
  }
  
  # Add UMAP coordinates if available
  if ("umap" %in% names(seurat_obj@reductions)) {
    umap_coords <- Embeddings(object = seurat_obj, reduction = "umap")
    expression_df$UMAP1 <- umap_coords[, 1]
    expression_df$UMAP2 <- umap_coords[, 2]
    expression_df$cluster <- Idents(object = seurat_obj)
  } else {
    warning("UMAP not present, generating dataframe without UMAP coordinates.\n")
  }
  
  return(expression_df)
}

# Joint Analyses
## joining and printing numbers
read_sample_files <- function(samples, directory, filtered = TRUE) {
  # Create an empty list to store data frames
  data_list <- list()
  
  all_data <- NULL
  
  for (sample in samples) {
    if (filtered) {
      # Define file types to look for when filtered = TRUE
      file_types <- c("SALVE", "GEX", "KLRB1")
      
      for (file_type in file_types) {
        file_path <- file.path(directory, paste0(sample, "_", file_type, "_filtered.csv"))
        
        # Check if file exists
        if (file.exists(file_path)) {
          # Read the file
          tryCatch({
            sample_data <- read.csv(file_path)
            
            sample_data$sample <- sample
            
            if (is.null(all_data)) {
              all_data <- sample_data
            } else {
              # Get common columns
              common_cols <- intersect(colnames(all_data), colnames(sample_data))
              all_data <- rbind(all_data[, common_cols], sample_data[, common_cols])
            }
            
            cat(sprintf("Successfully read %s\n", paste0(sample, "_", file_type, "_filtered.csv")))
          }, error = function(e) {
            cat(sprintf("Error reading %s: %s\n", file_path, conditionMessage(e)))
          })
        }
      }
    } else {
      # When filtered = FALSE, only look for one file type per sample
      file_path <- file.path(directory, paste0(sample, "_UMI_read_counts_full.csv"))
      
      # Check if file exists
      if (file.exists(file_path)) {
        # Read the file
        tryCatch({
          sample_data <- read.csv(file_path)
          
          sample_data$sample <- sample
          
          if (is.null(all_data)) {
            all_data <- sample_data
          } else {
            # Get common columns
            common_cols <- intersect(colnames(all_data), colnames(sample_data))
            all_data <- rbind(all_data[, common_cols], sample_data[, common_cols])
          }
          
          cat(sprintf("Successfully read %s\n", paste0(sample, "_UMI_read_counts_full.csv")))
        }, error = function(e) {
          cat(sprintf("Error reading %s: %s\n", file_path, conditionMessage(e)))
        })
      }
    }
  }
  
  if (!is.null(all_data)) {
    cat("Total cells loaded:", nrow(all_data), "\n")
    return(all_data)
  } else {
    warning("No data was loaded for any of the provided samples")
    return(NULL)
  }
}
## correlation analysis with plots of relationship between SALVE and SingleCell data
correlation_plots <- function(salve_data, single_cell_data, join = "inner") {
  if ("VISERcount" %in% colnames(salve_data)) {
    salve_data <- rename(salve_data, SALVEcount = VISERcount)
  }
  required_columns <- list(
    salve_data = c("SALVEcount", "cellID"),
    single_cell_data = c("SingleCellcount", "cellID")
  )
  
  lapply(names(required_columns), function(dataset) {
    missing_cols <- required_columns[[dataset]][!required_columns[[dataset]] %in% colnames(get(dataset))]
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns in %s: %s", dataset, paste(missing_cols, collapse = ", ")))
    }
  })
  
  join_func <- switch(join,
                      "inner" = inner_join,
                      "left" = left_join,
                      "right" = right_join,
                      "full" = full_join,
                      inner_join
  )

  jointTable <- join_func(salve_data, single_cell_data, by = "cellID") %>%
    select(-contains("UMAP")) %>%
    mutate(
      log1pSalveCount = log1p(SALVEcount),
      log2SalveCount = log2(SALVEcount),
      log1pSingleCell = log1p(SingleCellcount),
    )
  
  cat(
    "Correlation between 10X data and log1p SALVE data:\n",
    "  Pearson: ", round(cor(jointTable$log1pSingleCell, jointTable$log1pSalveCount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$log1pSingleCell, jointTable$log1pSalveCount, method="spearman"), 3),
    
    "\n\nCorrelation between 10X data and log2 SALVE data:\n",
    "  Pearson: ", round(cor(jointTable$log1pSingleCell, jointTable$log2SalveCount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$log1pSingleCell, jointTable$log2SalveCount, method="spearman"), 3), 
    
    "\n\nCorrelation between 10X counts and log2 SALVE data:\n",
    "  Pearson: ", round(cor(jointTable$SingleCellcount, jointTable$log2SalveCount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$SingleCellcount, jointTable$log2SalveCount, method="spearman"), 3),
    
    "\n\nCorrelation between 10X counts and SALVE counts:\n",
    "  Pearson: ", round(cor(jointTable$SingleCellcount, jointTable$SALVEcount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$SingleCellcount, jointTable$SALVEcount, method="spearman"), 3)
  )
  
  transformations <- list(
    SALVEcount = "Raw SALVE counts",
    log1pSalveCount = "Log1p SALVE",
    log2SalveCount = "Log2 SALVE"
  )
  
  plots <- lapply(names(transformations), function(trans) {
    ggplot(data = jointTable, aes_string(x = trans, y = "SingleCellcount")) +
      geom_point() +
      labs(
        x = transformations[[trans]],
        y = "SingleCell counts",
        title = paste("SingleCell vs", transformations[[trans]])
      ) +
      theme_minimal()
  })
  names(plots) <- names(transformations)
  
  return(list(
    data = jointTable,
    plots = plots
  ))
}
## compare host gene expression to viral expression
host_virus_correlation <- function(joint_salve_data, joint_rds, genes_list, viral_col,
                                   output.dir = "/projects/b1042/GoyalLab/egrody/extractedData/EGS018/joint/hostVirus/", 
                                   outputs = TRUE) {
  if (outputs && !dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  for (sample_name in names(joint_salve_data)) {
    current_df <- joint_salve_data[[sample_name]]
    seurat_obj <- joint_rds[[sample_name]]
    expression_df <- targetExpressionDF(seurat_obj, genes = genes_list) %>% select(-UMAP1, -UMAP2, -cluster)
    
    genes <- setdiff(colnames(expression_df), "cellID")
    if (!viral_col %in% colnames(current_df)) {
      cat("=== SAMPLE:", sample_name, "===\n")
      cat("Column", viral_col, "not found in data. Skipping.\n\n")
      next
    }
    salve <- current_df[[viral_col]]
    if (sd(salve) == 0) {
      cat("=== SAMPLE:", sample_name, "===\n")
      cat("SALVE has zero variance (all values identical). Skipping correlation analysis.\n\n")
      next
    }
    results <- data.frame()
    for (gene in genes) {
      gene_expr <- expression_df[[gene]]
      # Use cells where at least one is expressed
      cells_use <- (gene_expr > 0) | (salve > 0)
      n_cells <- sum(cells_use)
      if (n_cells < 10) {
        next
      }
      # Check for sufficient variance in both variables
      gene_subset <- gene_expr[cells_use]
      salve_subset <- salve[cells_use]
      
      if (sd(gene_subset) == 0 || sd(salve_subset) == 0) {
        next  # Skip genes with zero variance
      }
      
      # Spearman correlation
      cor_test <- cor.test(gene_expr[cells_use], 
                           salve[cells_use],
                           method = "spearman",
                           exact = FALSE)
      
      results <- rbind(results, data.frame(
        gene = gene,
        correlation = as.numeric(cor_test$estimate),
        p_value = cor_test$p.value,
        n_cells = n_cells,
        gene_positive = sum(gene_expr > 0),
        salve_positive = sum(salve > 0),
        both_positive = sum((gene_expr > 0) & (salve > 0))
      ))
    }
    # Adjust p-values for FDR
    results$p_adj <- p.adjust(results$p_value, method = "BH")
    # Sort by significance
    results <- results %>% 
      filter(p_adj < 0.05) %>%
      arrange(desc(abs(correlation)))
    
    # Print results
    cat("===", sample_name, "RESULTS (ordered by magnitude) ===\n\n")
    cat(sprintf("%-12s %10s %12s %12s %8s\n", 
                "Gene", "Corr", "P-value", "P-adj", "Sig"))
    cat(strrep("-", 65), "\n")
    
    for (i in 1:nrow(results)) {
      r <- results[i,]
      sig <- if (is.na(r$p_adj)) "NA" else if (r$p_adj < 0.001) "***" else if (r$p_adj < 0.01) "**" else if (r$p_adj < 0.05) "*" else "ns"
      cat(sprintf("%-12s %10.4f %12.2e %12.2e %8s\n",
                  r$gene, r$correlation, r$p_value, r$p_adj, sig))
    }
    cat(strrep("-", 65), "\n")
    
    # Save outputs if requested
    if (outputs) {
      # Add sample name column for context
      results$sample <- sample_name
      
      # Save CSV
      csv_file <- file.path(output.dir, paste0(sample_name, "_", viral_col, "_correlations.csv"))
      write.csv(results, csv_file, row.names = FALSE)
      
      # Create plots
      pdf_file <- file.path(output.dir, paste0(sample_name, "_", viral_col, "_correlation_barplots.pdf"))
      pdf(pdf_file, width = 10, height = 8)
      
      # Plot 1: Correlation coefficient barplot
      p1 <- ggplot(results, aes(x = reorder(gene, abs(correlation)), y = correlation, 
                                fill = correlation < 0)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                          labels = c("Positive", "Negative"),
                          name = "Direction") +
        labs(title = paste0(sample_name, " Gene-", viral_col, " Correlations"),
             subtitle = "Spearman correlation coefficients (FDR < 0.05)",
             x = "Gene",
             y = "Correlation Coefficient") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold", size = 14),
              axis.text = element_text(size = 10))
      print(p1)
      dev.off()
      
      pdf_file <- file.path(output.dir, paste0(sample_name, "_", viral_col, "_correlation_scatterplots.pdf"))
      pdf(pdf_file, width = 10, height = 8)
      # Plot 2: Scatter plots for top 6 genes
      top_genes <- head(results$gene, 6)
      
      for (gene in top_genes) {
        gene_expr <- expression_df[[gene]]
        plot_data <- data.frame(
          gene_expr = gene_expr,
          salve = salve
        )
        
        # Add jitter for better visualization of discrete values
        p2 <- ggplot(plot_data, aes(x = gene_expr, y = salve)) +
          geom_point(alpha = 0.3, size = 2) +
          geom_smooth(method = "lm", color = "#d73027", se = TRUE) +
          labs(title = paste(sample_name, "-", gene, "vs", viral_col),
               subtitle = sprintf("Spearman rho = %.3f, p-adj = %.2e",
                                  results$correlation[results$gene == gene],
                                  results$p_adj[results$gene == gene]),
               x = paste(gene, "expression"),
               y = viral_col) +
          theme_minimal() +
          theme(plot.title = element_text(face = "bold"))
        print(p2)
      }
      dev.off()
    }
  }
}

# Plotting
## plot a UMAP
plotUMAP <- function(data, colorby, title, output_dir, saveas, comparison = FALSE, color = "darkblue") {
  # Convert unquoted column name to string
  colorby_str <- deparse(substitute(colorby))
  
  # Create a copy of the data and arrange it by the color column (low to high)
  data_ordered <- data[order(data[[colorby_str]]), ]
  
  umap <- ggplot(data_ordered, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = .data[[colorby_str]]), size = 1, shape = 16) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.6)),
          legend.text = element_text(size = rel(0.6), angle = 30),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(title = title, color = "")
  
  if (comparison) {
    mid <- median(paint_umap$ratioSingleVISER)
    umap <- umap + scale_color_gradient2(midpoint = mid, low = "blue", mid = "gray93", high = "red")
  } else {
    umap <- umap + scale_color_gradient(low = "lightgrey", high = color)
  }
  ggsave(umap, file = paste0(output_dir, saveas),
         device = "pdf")
  
  return(umap)
}
## joint UMAP plot function
create_identity_umap <- function(data, title, output_file_name, output_dir = "/projects/b1042/GoyalLab/egrody/") {
  if ("log2ViserCount" %in% colnames(data)) {
    data <- rename(data, log1pSALVE = log2ViserCount)
  }
  # Check data is correct format
  required_columns <- c("log1pSALVE", "log1pSingleCell", "UMAP1", "UMAP2")
  missing_columns <- required_columns[!required_columns %in% colnames(data)]
  
  if (length(missing_columns) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_columns, collapse = ", ")))
  }
  
  # Create identity column
  paint_identities_umap <- data %>% 
    mutate(identity = case_when(
      log1pSALVE > 0 & log1pSingleCell > 0 ~ "both",
      log1pSALVE > 0 ~ "SALVE",
      log1pSingleCell > 0 ~ "SingleCell",
      log1pSALVE == 0 & log1pSingleCell == 0 ~ "neither"
    )) %>%
    mutate(identity = factor(identity))
  
  # Split data by identity
  neither_points <- paint_identities_umap %>% filter(identity == "neither")
  SALVE_points <- paint_identities_umap %>% filter(identity == "SALVE")
  both_points <- paint_identities_umap %>% filter(identity == "both")
  SingleCell_points <- paint_identities_umap %>% filter(identity == "SingleCell")
  
  identity_umap <- ggplot() +
    # Plot points in layers by identity
    geom_point(data = neither_points, 
               aes(x = UMAP1, y = UMAP2, color = identity),
               size = 1, shape = 16) +
    geom_point(data = SALVE_points, 
               aes(x = UMAP1, y = UMAP2, color = identity),
               size = 1, shape = 16) +
    geom_point(data = both_points, 
               aes(x = UMAP1, y = UMAP2, color = identity),
               size = 1, shape = 16) +
    geom_point(data = SingleCell_points, 
               aes(x = UMAP1, y = UMAP2, color = identity),
               size = 1, shape = 16) +
    scale_color_manual(values = c("both" = "purple", "neither" = "gray93", "SALVE" = "blue", "SingleCell" = "red")) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.6)),
          legend.text = element_text(size = rel(0.6), angle = 30),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(title = title, color = "Identity")
  
  ggsave(identity_umap, file = file.path(output_dir, output_file_name), width = 7, height = 7)
  
  return(identity_umap)
}
## barcode rank plot
barcodeRankPlot <- function(rawDataFolder, jointFullJoin, plotTitle = "Barcode Rank Plot", 
                            output_dir = "/projects/b1042/GoyalLab/egrody/", saveas = "barcodeRankPlot.svg") {
  # Check for either log1pSALVE or total_lessLTR columns
  has_log1p <- "log1pSALVE" %in% colnames(jointFullJoin)
  has_total <- "total_lessLTR" %in% colnames(jointFullJoin)
  if (!has_log1p && !has_total) {
    stop("Error: jointFullJoin must contain either a 'log1pSALVE' or 'total_lessLTR' column")
  }
  if (!"cellID" %in% colnames(jointFullJoin)) {
    stop("Error: jointFullJoin must contain a 'cellID' column")
  }
  
  # Determine which column to use
  if (has_log1p) {
    salve_col <- "log1pSALVE"
  } else {
    salve_col <- "total_lessLTR"
  }
  
  # Create a data frame with barcode and UMI count
  rawdata <- Read10X(rawDataFolder)
  umi_counts <- colSums(rawdata)
  plot_data <- data.frame(barcode = names(umi_counts), 
                          umi_count = umi_counts)
  plot_data <- plot_data %>% 
    arrange(desc(umi_count))
  
  # Filter for positive SALVE values
  # For log1pSALVE, use > 0; for total_lessLTR, use > 0 (raw counts)
  SALVE_cellids <- jointFullJoin %>% 
    filter(!!sym(salve_col) > 0) %>% 
    select(cellID)
  
  cat(sum(SALVE_cellids[[1]] %in% plot_data$barcode), "out of", nrow(SALVE_cellids), 
      "SALVEseq cellIDs are found in raw cellID list\n")
  
  plot_data$highlight <- plot_data$barcode %in% SALVE_cellids[[1]]
  
  # Count SALVE+ cellIDs with transcriptomes
  SALVE_in_rawdata <- plot_data %>% 
    filter(barcode %in% SALVE_cellids[[1]])
  
  num_full_transcriptomes <- sum(SALVE_in_rawdata$umi_count >= 100)
  cat(num_full_transcriptomes, "out of", nrow(SALVE_cellids), 
      "SALVEseq cellIDs have full transcriptomes (>= 100 UMI)\n")
  
  set.seed(123)
  plot_data <- plot_data %>%
    arrange(desc(umi_count)) %>%
    mutate(
      rank = row_number(),
      umi_count_adj = ifelse(umi_count <= 0, 0.1, umi_count)  # Replace 0 or negative with 0.1
    ) %>%
    slice_sample(n = min(20000, nrow(.))) # Update if necessary, ensure n doesn't exceed available rows
  
  # Create the plot
  umi_rank_plot <- ggplot(plot_data %>% arrange(highlight), aes(x = rank, y = umi_count_adj)) +
    geom_point(aes(color = highlight), alpha = 0.6) +
    scale_x_log10(
      labels = scales::label_number(),
      breaks = scales::breaks_log(n = 6)
    ) +
    scale_y_log10(
      labels = scales::label_number(),
      breaks = scales::breaks_log(n = 6)
    ) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    labs(x = "Rank of Barcode", 
         y = "UMI Count (Adjusted)",
         title = plotTitle) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Save the plot
  output_path <- file.path(output_dir, saveas)
  cat("Saving plot to:", output_path, "\n")
  ggsave(output_path, umi_rank_plot, width = 10, height = 8, device = "svg")
  
  return(umi_rank_plot)
}

# DEG Analysis
## perform comprehensive subset analysis
analyze_cell_subset <- function(seurat_obj, 
                                subset_cells_1, 
                                subset_cells_2 = NULL, 
                                min_pct = 0.1,
                                group_names = c("group1", "group2")) {
  
  # Input validation
  if (!is(seurat_obj, "Seurat")) {
    stop("First argument must be a Seurat object")
  }
  
  if (!is.vector(subset_cells_1)) {
    stop("subset_cells_1 must be a vector of cell IDs")
  }
  
  # Determine comparison mode
  two_group_mode <- !is.null(subset_cells_2)
  
  if (two_group_mode) {
    if (!is.vector(subset_cells_2)) {
      stop("subset_cells_2 must be a vector of cell IDs")
    }
    cat("Running two-group comparison mode\n")
  } else {
    cat("Running single-group comparison mode (subset vs rest)\n")
  }
  
  # Convert cell names to character vectors
  subset_cells_1 <- as.character(subset_cells_1)
  if (two_group_mode) {
    subset_cells_2 <- as.character(subset_cells_2)
    
    # Check for overlap between groups
    overlap <- intersect(subset_cells_1, subset_cells_2)
    if (length(overlap) > 0) {
      stop(sprintf("Found %d cells present in both subset_cells_1 and subset_cells_2. Groups must be mutually exclusive.",
                   length(overlap)))
    }
  }
  
  # Validate subset_cells_1
  valid_cells_1 <- subset_cells_1 %in% colnames(seurat_obj)
  if (!any(valid_cells_1)) {
    stop("None of the provided cell IDs in subset_cells_1 were found in the Seurat object")
  }
  
  if (!all(valid_cells_1)) {
    warning(sprintf("%d cells from subset_cells_1 were not found in the Seurat object",
                    sum(!valid_cells_1)))
    subset_cells_1 <- subset_cells_1[valid_cells_1]
  }
  
  # Validate subset_cells_2 if provided
  if (two_group_mode) {
    valid_cells_2 <- subset_cells_2 %in% colnames(seurat_obj)
    if (!any(valid_cells_2)) {
      stop("None of the provided cell IDs in subset_cells_2 were found in the Seurat object")
    }
    
    if (!all(valid_cells_2)) {
      warning(sprintf("%d cells from subset_cells_2 were not found in the Seurat object",
                      sum(!valid_cells_2)))
      subset_cells_2 <- subset_cells_2[valid_cells_2]
    }
  }
  
  # Print diagnostic information
  cat(sprintf("Found %d cells in %s\n", length(subset_cells_1), group_names[1]))
  if (two_group_mode) {
    cat(sprintf("Found %d cells in %s\n", length(subset_cells_2), group_names[2]))
  } else {
    n_rest <- ncol(seurat_obj) - length(subset_cells_1)
    cat(sprintf("Comparing against %d cells in rest of dataset\n", n_rest))
  }
  
  # Add metadata column for group membership
  if (two_group_mode) {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1],
      ifelse(colnames(seurat_obj) %in% subset_cells_2, group_names[2], "excluded")
    )
    # Subset to only include cells in one of the two groups
    seurat_obj <- subset(seurat_obj, cells = c(subset_cells_1, subset_cells_2))
    cat(sprintf("Analysis will include %d total cells\n", ncol(seurat_obj)))
  } else {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1], group_names[2]
    )
  }
  
  # Check for cluster information
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    stop("No 'seurat_clusters' column found in metadata. Please run FindClusters() first")
  }
  
  # Get unique clusters for each group
  clusters_group1 <- unique(seurat_obj$seurat_clusters[seurat_obj$comparison_group == group_names[1]])
  clusters_group2 <- unique(seurat_obj$seurat_clusters[seurat_obj$comparison_group == group_names[2]])
  
  cat(sprintf("%s cells are found in %d clusters: %s\n",
              group_names[1], length(clusters_group1),
              paste(sort(clusters_group1), collapse = ", ")))
  cat(sprintf("%s cells are found in %d clusters: %s\n",
              group_names[2], length(clusters_group2),
              paste(sort(clusters_group2), collapse = ", ")))
  
  # Get overlapping clusters
  overlap_clusters <- intersect(clusters_group1, clusters_group2)
  cat(sprintf("Overlapping clusters: %d (%s)\n",
              length(overlap_clusters),
              ifelse(length(overlap_clusters) > 0, 
                     paste(sort(overlap_clusters), collapse = ", "),
                     "none")))
  
  # Initialize results list
  results <- list()
  results$metadata <- list(
    mode = ifelse(two_group_mode, "two_group", "single_group"),
    comparison_summary = sprintf("Comparing %s (n=%d) vs %s (n=%d)", 
                                 group_names[1], length(subset_cells_1),
                                 group_names[2], 
                                 ifelse(two_group_mode, length(subset_cells_2), 
                                        ncol(seurat_obj) - length(subset_cells_1))),
    group1_name = group_names[1],
    group2_name = group_names[2],
    n_cells_group1 = length(subset_cells_1),
    n_cells_group2 = ifelse(two_group_mode, length(subset_cells_2), 
                            ncol(seurat_obj) - length(subset_cells_1)),
    group_names = group_names,
    overlap_clusters = overlap_clusters
  )
  
  # 1. Overall differential expression (across all clusters)
  cat("\n=== Overall Differential Expression Analysis ===\n")
  Idents(seurat_obj) <- factor(seurat_obj$comparison_group, levels = group_names)
  
  results$overall_de <- tryCatch({
    de_results <- FindMarkers(seurat_obj,
                              ident.1 = group_names[1],
                              ident.2 = group_names[2],
                              min.pct = min_pct,
                              test.use = "wilcox",
                              verbose = FALSE)
    
    cat(sprintf("Found %d differentially expressed genes (all p-values)\n", nrow(de_results)))
    de_results
  }, error = function(e) {
    warning(sprintf("Error in overall DE analysis: %s", e$message))
    NULL
  })
  
  # Filter for significant overall markers
  results$overall_de_significant <- tryCatch({
    if (!is.null(results$overall_de) && nrow(results$overall_de) > 0) {
      sig_results <- results$overall_de[
        abs(results$overall_de$avg_log2FC) > 0.5 & 
          results$overall_de$p_val_adj < 0.01,
      ]
      
      if (nrow(sig_results) > 0) {
        n_up <- sum(sig_results$avg_log2FC > 0)
        n_down <- nrow(sig_results) - n_up
        cat(sprintf("Significant markers: %d total (%d up in %s, %d down)\n",
                    nrow(sig_results), n_up, group_names[1], n_down))
        sig_results$gene <- rownames(sig_results)
        sig_results
      } else {
        cat("No markers passed significance thresholds\n")
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error filtering overall markers: ", e$message)
    NULL
  })
  
  # 2. Cluster-specific analysis (only for overlapping clusters)
  if (length(overlap_clusters) > 0) {
    cat("\n=== Cluster-Specific Differential Expression Analysis ===\n")
    cluster_results <- list()
    
    for(cluster in overlap_clusters) {
      # Get cells from this cluster for each group
      cluster_cells <- seurat_obj$seurat_clusters == cluster
      group1_in_cluster <- sum(seurat_obj$comparison_group[cluster_cells] == group_names[1])
      group2_in_cluster <- sum(seurat_obj$comparison_group[cluster_cells] == group_names[2])
      
      # Require at least 3 cells per group for statistical power
      if (group1_in_cluster < 3 || group2_in_cluster < 3) {
        warning(sprintf("Skipping cluster %s: insufficient cells (%s: %d, %s: %d)\n",
                        cluster, group_names[1], group1_in_cluster, 
                        group_names[2], group2_in_cluster))
        next
      }
      
      # Create a temporary Seurat object for this cluster
      cluster_obj <- subset(seurat_obj, seurat_clusters == cluster)
      Idents(cluster_obj) <- factor(cluster_obj$comparison_group, levels = group_names)
      
      # Perform DE analysis within this cluster
      tryCatch({
        de_results <- FindMarkers(cluster_obj,
                                  ident.1 = group_names[1],
                                  ident.2 = group_names[2],
                                  min.pct = min_pct,
                                  test.use = "wilcox",
                                  verbose = FALSE)
        
        cluster_results[[paste0("cluster_", cluster)]] <- de_results
        
        cat(sprintf("Cluster %s: %d DEGs (%s: %d cells, %s: %d cells)\n",
                    cluster, nrow(de_results), 
                    group_names[1], group1_in_cluster,
                    group_names[2], group2_in_cluster))
        
      }, error = function(e) {
        warning(sprintf("Error in DE analysis for cluster %s: %s", cluster, e$message))
      })
    }
    
    results$cluster_specific <- cluster_results
    
    # Filter cluster-specific results for significance
    results$cluster_specific_significant <- tryCatch({
      if (length(cluster_results) > 0) {
        significant_by_cluster <- lapply(names(cluster_results), function(cluster_name) {
          de_results <- cluster_results[[cluster_name]]
          
          if (!is.null(de_results) && nrow(de_results) > 0) {
            sig_results <- de_results[
              abs(de_results$avg_log2FC) > 0.5 & 
                de_results$p_val_adj < 0.01,
            ]
            
            if (nrow(sig_results) > 0) {
              sig_results$gene <- rownames(sig_results)
              sig_results
            } else {
              NULL
            }
          } else {
            NULL
          }
        })
        
        names(significant_by_cluster) <- names(cluster_results)
        significant_by_cluster <- significant_by_cluster[!sapply(significant_by_cluster, is.null)]
        
        if (length(significant_by_cluster) > 0) {
          total_sig <- sum(sapply(significant_by_cluster, nrow))
          cat(sprintf("\nCluster-specific significant markers: %d genes across %d clusters\n",
                      total_sig, length(significant_by_cluster)))
          
          for (cluster_name in names(significant_by_cluster)) {
            n_markers <- nrow(significant_by_cluster[[cluster_name]])
            n_up <- sum(significant_by_cluster[[cluster_name]]$avg_log2FC > 0)
            n_down <- n_markers - n_up
            cat(sprintf("  %s: %d markers (%d up in %s, %d down)\n", 
                        cluster_name, n_markers, group_names[1], n_up, n_down))
          }
          
          significant_by_cluster
        } else {
          cat("\nNo significant cluster-specific markers found\n")
          NULL
        }
      } else {
        NULL
      }
    }, error = function(e) {
      warning("Error filtering cluster-specific markers: ", e$message)
      NULL
    })
  } else {
    cat("\n=== Cluster-Specific Analysis Skipped ===\n")
    cat("No overlapping clusters between groups\n")
    results$cluster_specific <- NULL
    results$cluster_specific_significant <- NULL
  }
  
  # 3. Pseudo-bulk expression profiles
  cat("\n=== Pseudo-bulk Expression ===\n")
  results$pseudobulk_expression <- tryCatch({
    pb_expr <- AggregateExpression(seurat_obj,
                                   group.by = "comparison_group",
                                   assays = "RNA",
                                   return.seurat = FALSE)$RNA
    cat(sprintf("Calculated pseudo-bulk expression for %d genes\n", nrow(pb_expr)))
    pb_expr
  }, error = function(e) {
    warning("Error calculating pseudo-bulk expression: ", e$message)
    NULL
  })
  
  # Print final summary
  cat("\n=== Analysis Summary ===\n")
  cat(sprintf("Comparison: %s\n", results$metadata$comparison_summary))
  cat(sprintf("Mode: %s\n", results$metadata$mode))
  
  return(results)
}
## analyze_cell_subset, without using clusters
analyze_cell_subset_unclustered <- function(seurat_obj, 
                                            subset_cells_1, 
                                            subset_cells_2 = NULL,
                                            min_pct = 0.1,
                                            group_names = c("group1", "group2")) {
  # Check package dependencies
  rlang::check_installed("Seurat", reason = "to perform single-cell analysis")
  rlang::check_installed("dplyr", reason = "for data manipulation")
  
  require(Seurat)
  require(dplyr)
  
  # Input validation
  if (!is(seurat_obj, "Seurat")) {
    stop("First argument must be a Seurat object")
  }
  
  if (!is.vector(subset_cells_1)) {
    stop("subset_cells_1 must be a vector of cell IDs")
  }
  
  # Determine comparison mode
  two_group_mode <- !is.null(subset_cells_2)
  
  if (two_group_mode) {
    if (!is.vector(subset_cells_2)) {
      stop("subset_cells_2 must be a vector of cell IDs")
    }
    cat("Running two-group comparison mode (unclustered)\n")
  } else {
    cat("Running single-group comparison mode (subset vs rest, unclustered)\n")
  }
  
  # Convert cell names to character vectors
  subset_cells_1 <- as.character(subset_cells_1)
  if (two_group_mode) {
    subset_cells_2 <- as.character(subset_cells_2)
    
    # Check for overlap between groups
    overlap <- intersect(subset_cells_1, subset_cells_2)
    if (length(overlap) > 0) {
      stop(sprintf("Found %d cells present in both subset_cells_1 and subset_cells_2. Groups must be mutually exclusive.",
                   length(overlap)))
    }
  }
  
  # Validate subset_cells_1
  valid_cells_1 <- subset_cells_1 %in% colnames(seurat_obj)
  if (!any(valid_cells_1)) {
    stop("None of the provided cell IDs in subset_cells_1 were found in the Seurat object")
  }
  
  if (!all(valid_cells_1)) {
    n_invalid <- sum(!valid_cells_1)
    n_valid <- sum(valid_cells_1)
    warning(sprintf("%d cells from subset_cells_1 were not found in the Seurat object. %d valid cells remain.",
                    n_invalid, n_valid))
    
    # Show some examples of invalid cells for debugging
    if (n_invalid <= 5) {
      cat(sprintf("Invalid cell IDs: %s\n", paste(subset_cells_1[!valid_cells_1], collapse = ", ")))
    }
    
    subset_cells_1 <- subset_cells_1[valid_cells_1]
  }
  
  # Early check for insufficient cells
  if (length(subset_cells_1) < 3) {
    stop(sprintf("Insufficient valid cells in subset_cells_1: found %d, need at least 3 for statistical analysis.\nConsider checking your cell ID format or filtering criteria.",
                 length(subset_cells_1)))
  }
  
  # Validate subset_cells_2 if provided
  if (two_group_mode) {
    valid_cells_2 <- subset_cells_2 %in% colnames(seurat_obj)
    if (!any(valid_cells_2)) {
      stop("None of the provided cell IDs in subset_cells_2 were found in the Seurat object")
    }
    
    if (!all(valid_cells_2)) {
      n_invalid <- sum(!valid_cells_2)
      n_valid <- sum(valid_cells_2)
      warning(sprintf("%d cells from subset_cells_2 were not found in the Seurat object. %d valid cells remain.",
                      n_invalid, n_valid))
      
      # Show some examples of invalid cells for debugging
      if (n_invalid <= 5) {
        cat(sprintf("Invalid cell IDs: %s\n", paste(subset_cells_2[!valid_cells_2], collapse = ", ")))
      }
      
      subset_cells_2 <- subset_cells_2[valid_cells_2]
    }
    
    # Early check for insufficient cells
    if (length(subset_cells_2) < 3) {
      stop(sprintf("Insufficient valid cells in subset_cells_2: found %d, need at least 3 for statistical analysis.\nConsider checking your cell ID format or filtering criteria.",
                   length(subset_cells_2)))
    }
  }
  
  # Print diagnostic information
  cat(sprintf("\nDataset information:\n"))
  cat(sprintf("  - Total cells in dataset: %d\n", ncol(seurat_obj)))
  cat(sprintf("  - %s: %d cells (%.1f%% of total)\n", 
              group_names[1], 
              length(subset_cells_1),
              100 * length(subset_cells_1) / ncol(seurat_obj)))
  
  if (two_group_mode) {
    cat(sprintf("  - %s: %d cells (%.1f%% of total)\n",
                group_names[2],
                length(subset_cells_2),
                100 * length(subset_cells_2) / ncol(seurat_obj)))
  } else {
    n_rest <- ncol(seurat_obj) - length(subset_cells_1)
    cat(sprintf("  - %s: %d cells (%.1f%% of total)\n",
                group_names[2],
                n_rest,
                100 * n_rest / ncol(seurat_obj)))
  }
  
  # Add metadata column for group membership
  if (two_group_mode) {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1],
      ifelse(colnames(seurat_obj) %in% subset_cells_2, group_names[2], "excluded")
    )
    # Subset to only include cells in one of the two groups
    seurat_obj <- subset(seurat_obj, cells = c(subset_cells_1, subset_cells_2))
    cat(sprintf("  - Analysis will include %d total cells\n", ncol(seurat_obj)))
  } else {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1], group_names[2]
    )
  }
  
  # Count cells in each group
  n_group1 <- sum(seurat_obj$comparison_group == group_names[1])
  n_group2 <- sum(seurat_obj$comparison_group == group_names[2])
  
  cat(sprintf("\nComparing:\n"))
  cat(sprintf("  - %s: %d cells\n", group_names[1], n_group1))
  cat(sprintf("  - %s: %d cells\n", group_names[2], n_group2))
  
  # Initialize results list
  results <- list()
  results$metadata <- list(
    mode = ifelse(two_group_mode, "two_group", "single_group"),
    n_cells_group1 = n_group1,
    n_cells_group2 = n_group2,
    group_names = group_names,
    analysis_type = "unclustered"
  )
  
  # Check if we have enough cells for comparison
  if (n_group1 < 3 || n_group2 < 3) {
    stop(sprintf("Insufficient cells for comparison (%s: %d, %s: %d). Need at least 3 in each group.",
                 group_names[1], n_group1, group_names[2], n_group2))
  }
  
  # 1. Perform unclustered differential expression analysis
  cat("\n=== Differential Expression Analysis ===\n")
  
  # Set identities for comparison
  Idents(seurat_obj) <- factor(seurat_obj$comparison_group, levels = group_names)
  
  # Perform differential expression
  results$de_markers <- tryCatch({
    de_results <- FindMarkers(seurat_obj,
                              ident.1 = group_names[1],
                              ident.2 = group_names[2],
                              min.pct = min_pct,
                              test.use = "wilcox",
                              verbose = FALSE)
    
    # Add gene column
    de_results$gene <- rownames(de_results)
    
    cat(sprintf("Successfully tested %d genes\n", nrow(de_results)))
    cat(sprintf("  - Significant at p_val_adj < 0.05: %d genes\n", 
                sum(de_results$p_val_adj < 0.05)))
    cat(sprintf("  - Significant at p_val_adj < 0.01: %d genes\n",
                sum(de_results$p_val_adj < 0.01)))
    
    de_results
    
  }, error = function(e) {
    stop(sprintf("Error in differential expression analysis: %s", e$message))
  })
  
  # 2. Filter for significant markers
  cat("\nFiltering for significant markers (|avg_log2FC| > 0.5, p_val_adj < 0.01)...\n")
  
  results$de_markers_significant <- tryCatch({
    if (!is.null(results$de_markers) && nrow(results$de_markers) > 0) {
      sig_markers <- results$de_markers[
        abs(results$de_markers$avg_log2FC) > 0.5 & 
          results$de_markers$p_val_adj < 0.01,
      ]
      
      if (nrow(sig_markers) > 0) {
        # Count up and down
        n_up <- sum(sig_markers$avg_log2FC > 0)
        n_down <- sum(sig_markers$avg_log2FC < 0)
        
        cat(sprintf("Found %d significant markers\n", nrow(sig_markers)))
        cat(sprintf("  - Upregulated in %s: %d genes\n", group_names[1], n_up))
        cat(sprintf("  - Downregulated in %s: %d genes\n", group_names[1], n_down))
        
        sig_markers
      } else {
        cat("No markers passed significance thresholds\n")
        NULL
      }
    } else {
      cat("No differential expression results to filter\n")
      NULL
    }
  }, error = function(e) {
    warning("Error filtering significant markers: ", e$message)
    NULL
  })
  
  # 3. Additional filtering: moderate fold change threshold
  results$de_markers_moderate <- tryCatch({
    if (!is.null(results$de_markers) && nrow(results$de_markers) > 0) {
      moderate_markers <- results$de_markers[
        abs(results$de_markers$avg_log2FC) > 0.25 & 
          results$de_markers$p_val_adj < 0.05,
      ]
      
      if (nrow(moderate_markers) > 0) {
        cat(sprintf("\nAlternative threshold (|avg_log2FC| > 0.25, p_val_adj < 0.05):\n"))
        cat(sprintf("  - Found %d markers\n", nrow(moderate_markers)))
        moderate_markers
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error applying moderate thresholds: ", e$message)
    NULL
  })
  
  # 4. Calculate pseudo-bulk expression profiles
  cat("\n=== Pseudo-bulk Expression ===\n")
  results$pseudobulk_expression <- tryCatch({
    pb_expr <- AggregateExpression(seurat_obj,
                                   group.by = "comparison_group",
                                   assays = "RNA",
                                   return.seurat = FALSE)$RNA
    cat(sprintf("Calculated pseudo-bulk expression for %d genes\n", nrow(pb_expr)))
    pb_expr
  }, error = function(e) {
    warning("Error calculating pseudo-bulk expression: ", e$message)
    NULL
  })
  
  # Summary of analysis
  cat("\n=== Analysis Summary ===\n")
  cat(sprintf("Mode: %s\n", results$metadata$mode))
  cat(sprintf("Analysis type: %s\n", results$metadata$analysis_type))
  cat(sprintf("Groups: %s (%d cells) vs %s (%d cells)\n",
              group_names[1], n_group1,
              group_names[2], n_group2))
  cat(sprintf("Total genes tested: %d\n", 
              ifelse(!is.null(results$de_markers), 
                     nrow(results$de_markers), 0)))
  cat(sprintf("Significant markers (strict): %s\n",
              ifelse(!is.null(results$de_markers_significant),
                     sprintf("%d genes", nrow(results$de_markers_significant)),
                     "None")))
  cat(sprintf("Significant markers (moderate): %s\n",
              ifelse(!is.null(results$de_markers_moderate),
                     sprintf("%d genes", nrow(results$de_markers_moderate)),
                     "None")))
  
  if (!is.null(results$de_markers_significant) && 
      nrow(results$de_markers_significant) > 0) {
    sig_markers <- results$de_markers_significant
    n_up <- sum(sig_markers$avg_log2FC > 0)
    n_down <- sum(sig_markers$avg_log2FC < 0)
    
    cat(sprintf("  - Upregulated in %s: %d\n", group_names[1], n_up))
    cat(sprintf("  - Downregulated in %s: %d\n", group_names[1], n_down))
    
    # Report top markers only if they exist
    if (n_up > 0) {
      top_up <- head(sig_markers[order(-sig_markers$avg_log2FC), ], 5)
      cat(sprintf("\nTop markers upregulated in %s:\n", group_names[1]))
      for (i in 1:nrow(top_up)) {
        cat(sprintf("  %d. %s (log2FC = %.2f, p_adj = %.2e)\n",
                    i, top_up$gene[i], top_up$avg_log2FC[i], top_up$p_val_adj[i]))
      }
    }
    
    if (n_down > 0) {
      top_down <- head(sig_markers[order(sig_markers$avg_log2FC), ], 5)
      cat(sprintf("\nTop markers downregulated in %s:\n", group_names[1]))
      for (i in 1:nrow(top_down)) {
        cat(sprintf("  %d. %s (log2FC = %.2f, p_adj = %.2e)\n",
                    i, top_down$gene[i], top_down$avg_log2FC[i], top_down$p_val_adj[i]))
      }
    }
  }
  
  cat("\n")
  
  return(results)
}
## visualize subset analysis results
plot_subset_analysis <- function(seurat_obj, analysis_results, subset_cells, save = TRUE, output_dir = "/projects/b1042/GoyalLab/egrody/") {
  require(ggplot2)
  
  if (!is(seurat_obj, "Seurat")) {
    stop("First argument must be a Seurat object")
  }
  
  if (!is.list(analysis_results)) {
    stop("analysis_results must be the output from analyze_cell_subset()")
  }
  
  plots <- list()
  
  # 1. UMAP highlighting subset cells
  if ("umap" %in% names(seurat_obj@reductions)) {
    seurat_obj$in_subset <- colnames(seurat_obj) %in% subset_cells
    plots$umap <- DimPlot(seurat_obj, 
                          group.by = "in_subset", 
                          label = TRUE) +
      ggtitle("Subset Cells in UMAP")
  } else {
    warning("No UMAP reduction found. Skipping UMAP plot.")
  }
  
  # 2. Expression heatmap of top markers
  if (!is.null(analysis_results$conserved_markers) && 
      nrow(analysis_results$conserved_markers) > 0) {
    top_conserved <- head(rownames(analysis_results$conserved_markers), 20)
    plots$heatmap <- DoHeatmap(seurat_obj, 
                               features = top_conserved,
                               group.by = c("in_subset", "seurat_clusters"))
  }
  
  # 3. Violin plots for top markers by cluster
  if (!is.null(analysis_results$top_genes)) {
    for(cluster in names(analysis_results$top_genes)) {
      top_genes <- analysis_results$top_genes[[cluster]]
      if (length(top_genes) > 0) {
        plots[[paste0("violin_", cluster)]] <- VlnPlot(seurat_obj,
                                                       features = head(top_genes, 5),
                                                       split.by = "in_subset",
                                                       group.by = "seurat_clusters")
      }
    }
  }
  
  if (save == TRUE) {
    # Create output directory if it doesn't exist
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Function to safely save a plot
    save_plot <- function(plot, filename, w = width, h = height) {
      if (!is.null(plot)) {
        ggsave(
          filename = file.path(output_dir, filename),
          plot = plot,
          width = w,
          height = h
        )
        return(TRUE)
      }
      return(FALSE)
    }
    
    # Initialize counters for summary
    saved_plots <- character()
    
    # 1. Save UMAP plot
    if (!is.null(plots$umap)) {
      if (save_plot(plots$umap, "umap_subset_cells.png", w = 8, h = 8)) {
        saved_plots <- c(saved_plots, "UMAP plot")
      }
    }
    
    # 2. Save heatmap
    if (!is.null(plots$heatmap)) {
      if (save_plot(plots$heatmap, "marker_heatmap.png", w = 12, h = 10)) {
        saved_plots <- c(saved_plots, "Marker heatmap")
      }
    }
    
    # 3. Save violin plots
    #violin_plots <- names(plots)[grep("^violin_", names(plots))]
    #if (length(violin_plots) > 0) {
    #  # Create a directory for violin plots
    #  violin_dir <- file.path(output_dir, "violin_plots")
    #  dir.create(violin_dir, showWarnings = FALSE)
    #  
    #  for (plot_name in violin_plots) {
    #    cluster <- gsub("violin_", "", plot_name)
    #    filename <- file.path("violin_plots", sprintf("cluster_%s_markers.png", cluster))
    #    if (save_plot(plots[[plot_name]], filename, w = 12, h = 8)) {
    #      saved_plots <- c(saved_plots, sprintf("Violin plot for cluster %s", cluster))
    #    }
    #  }
    print("Plots saved")
    }
  
  return(plots)
}
## save analysis results to files
save_subset_analysis <- function(analysis_results, output_dir) {
  require(openxlsx)
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 0. Save overall differential expression results
  if (!is.null(analysis_results$overall_de) && nrow(analysis_results$overall_de) > 0) {
    wb <- createWorkbook()
    
    # All results
    de_results <- analysis_results$overall_de
    de_results$gene <- rownames(de_results)
    addWorksheet(wb, "all_genes")
    writeData(wb, "all_genes", de_results)
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "overall_DE.xlsx"), overwrite = TRUE)
  }
  
  # 0b. Save overall SIGNIFICANT markers
  if (!is.null(analysis_results$overall_de_significant) && 
      nrow(analysis_results$overall_de_significant) > 0) {
    
    sig_results <- analysis_results$overall_de_significant
    
    # Split into upregulated and downregulated
    upregulated <- sig_results[sig_results$avg_log2FC > 0, ]
    downregulated <- sig_results[sig_results$avg_log2FC < 0, ]
    
    # Sort by fold change
    if (nrow(upregulated) > 0) {
      upregulated <- upregulated[order(-upregulated$avg_log2FC), ]
    }
    if (nrow(downregulated) > 0) {
      downregulated <- downregulated[order(downregulated$avg_log2FC), ]
    }
    
    wb <- createWorkbook()
    
    # Add upregulated genes tab
    if (nrow(upregulated) > 0) {
      addWorksheet(wb, "upregulated")
      writeData(wb, "upregulated", upregulated)
    }
    
    # Add downregulated genes tab
    if (nrow(downregulated) > 0) {
      addWorksheet(wb, "downregulated")
      writeData(wb, "downregulated", downregulated)
    }
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "overall_significant.xlsx"), overwrite = TRUE)
    
    # Print summary
    cat(sprintf("\nOverall significant markers saved:\n"))
    cat(sprintf("  - Upregulated genes: %d\n", nrow(upregulated)))
    cat(sprintf("  - Downregulated genes: %d\n", nrow(downregulated)))
  }
  
  # 1. Save cluster-specific differential expression results
  if (length(analysis_results$cluster_specific) > 0) {
    # Create a workbook for cluster-specific results
    wb <- createWorkbook()
    
    for (cluster_name in names(analysis_results$cluster_specific)) {
      de_results <- analysis_results$cluster_specific[[cluster_name]]
      if (!is.null(de_results) && nrow(de_results) > 0) {
        # Add cluster results as a worksheet
        addWorksheet(wb, cluster_name)
        # Add data with gene names as a column
        de_results$gene <- rownames(de_results)
        writeData(wb, cluster_name, de_results)
      }
    }
    
    # Save the workbook
    saveWorkbook(wb, file.path(output_dir, "cluster_specific_DE.xlsx"), overwrite = TRUE)
  }
  
  # 1b. Save cluster-specific SIGNIFICANT markers
  if (!is.null(analysis_results$cluster_specific_significant) &&
      length(analysis_results$cluster_specific_significant) > 0) {
    
    wb <- createWorkbook()
    
    for (cluster_name in names(analysis_results$cluster_specific_significant)) {
      sig_results <- analysis_results$cluster_specific_significant[[cluster_name]]
      
      if (!is.null(sig_results) && nrow(sig_results) > 0) {
        # Split into upregulated and downregulated
        upregulated <- sig_results[sig_results$avg_log2FC > 0, ]
        downregulated <- sig_results[sig_results$avg_log2FC < 0, ]
        
        # Sort by fold change
        if (nrow(upregulated) > 0) {
          upregulated <- upregulated[order(-upregulated$avg_log2FC), ]
        }
        if (nrow(downregulated) > 0) {
          downregulated <- downregulated[order(downregulated$avg_log2FC), ]
        }
        
        # Add worksheet for upregulated genes
        if (nrow(upregulated) > 0) {
          sheet_name_up <- paste0(cluster_name, "_up")
          # Truncate sheet name if too long (Excel limit is 31 characters)
          if (nchar(sheet_name_up) > 31) {
            sheet_name_up <- substr(sheet_name_up, 1, 31)
          }
          addWorksheet(wb, sheet_name_up)
          writeData(wb, sheet_name_up, upregulated)
        }
        
        # Add worksheet for downregulated genes
        if (nrow(downregulated) > 0) {
          sheet_name_down <- paste0(cluster_name, "_down")
          # Truncate sheet name if too long
          if (nchar(sheet_name_down) > 31) {
            sheet_name_down <- substr(sheet_name_down, 1, 31)
          }
          addWorksheet(wb, sheet_name_down)
          writeData(wb, sheet_name_down, downregulated)
        }
      }
    }
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "cluster_specific_significant.xlsx"), overwrite = TRUE)
  }
  
  # 2. Save conserved markers
  if (!is.null(analysis_results$conserved_markers) && 
      nrow(analysis_results$conserved_markers) > 0) {
    # Add gene names as a column
    conserved_markers <- analysis_results$conserved_markers
    conserved_markers$gene <- rownames(conserved_markers)
    
    # Create a new workbook for conserved markers
    wb <- createWorkbook()
    addWorksheet(wb, "conserved_markers")
    writeData(wb, "conserved_markers", conserved_markers)
    saveWorkbook(wb, file.path(output_dir, "conserved_markers.xlsx"), overwrite = TRUE)
  }
  
  # 2b. Save significant conserved markers
  if (!is.null(analysis_results$conserved_markers_significant) && 
      nrow(analysis_results$conserved_markers_significant) > 0) {
    # Add gene names as a column
    conserved_markers_sig <- analysis_results$conserved_markers_significant
    conserved_markers_sig$gene <- rownames(conserved_markers_sig)
    
    # Determine which fold change column to use for sorting
    fc_column <- if ("mean_log2FC" %in% colnames(conserved_markers_sig)) {
      "mean_log2FC"
    } else if ("avg_log2FC" %in% colnames(conserved_markers_sig)) {
      "avg_log2FC"
    } else {
      stop("No log2FC column found in conserved_markers_significant")
    }
    
    # Split into upregulated and downregulated
    upregulated <- conserved_markers_sig[conserved_markers_sig[[fc_column]] > 0, ]
    downregulated <- conserved_markers_sig[conserved_markers_sig[[fc_column]] < 0, ]
    
    # Sort by fold change (descending for up, ascending for down)
    if (nrow(upregulated) > 0) {
      upregulated <- upregulated[order(-upregulated[[fc_column]]), ]
    }
    if (nrow(downregulated) > 0) {
      downregulated <- downregulated[order(downregulated[[fc_column]]), ]
    }
    
    # Create a new workbook for significant conserved markers
    wb <- createWorkbook()
    
    # Add upregulated genes tab
    if (nrow(upregulated) > 0) {
      addWorksheet(wb, "upregulated")
      writeData(wb, "upregulated", upregulated)
    }
    
    # Add downregulated genes tab
    if (nrow(downregulated) > 0) {
      addWorksheet(wb, "downregulated")
      writeData(wb, "downregulated", downregulated)
    }
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "conserved_markers_significant.xlsx"), overwrite = TRUE)
    
    # Print summary
    cat(sprintf("\nSignificant conserved markers saved:\n"))
    cat(sprintf("  - Upregulated genes: %d\n", nrow(upregulated)))
    cat(sprintf("  - Downregulated genes: %d\n", nrow(downregulated)))
  }
  
  # 3. Save average expression profiles
  if (!is.null(analysis_results$avg_expression)) {
    avg_exp <- analysis_results$avg_expression
    avg_exp_df <- as.data.frame(avg_exp)
    avg_exp_df$gene <- rownames(avg_exp_df)
    
    wb <- createWorkbook()
    addWorksheet(wb, "average_expression")
    writeData(wb, "average_expression", avg_exp_df)
    saveWorkbook(wb, file.path(output_dir, "average_expression.xlsx"), overwrite = TRUE)
  }
  
  # 4. Save pseudobulk expression profiles if present
  if (!is.null(analysis_results$pseudobulk_expression)) {
    pseudobulk_exp <- analysis_results$pseudobulk_expression
    pseudobulk_exp_df <- as.data.frame(pseudobulk_exp)
    pseudobulk_exp_df$gene <- rownames(pseudobulk_exp_df)
    
    wb <- createWorkbook()
    addWorksheet(wb, "pseudobulk_expression")
    writeData(wb, "pseudobulk_expression", pseudobulk_exp_df)
    saveWorkbook(wb, file.path(output_dir, "pseudobulk_expression.xlsx"), overwrite = TRUE)
  }
  
  # 5. Save top genes list
  if (length(analysis_results$top_genes) > 0) {
    # Convert list to data frame
    max_genes <- max(sapply(analysis_results$top_genes, length))
    top_genes_df <- data.frame(matrix(NA, nrow = max_genes, 
                                      ncol = length(analysis_results$top_genes)))
    colnames(top_genes_df) <- names(analysis_results$top_genes)
    
    for (cluster in names(analysis_results$top_genes)) {
      top_genes_df[[cluster]] <- c(analysis_results$top_genes[[cluster]],
                                   rep(NA, max_genes - length(analysis_results$top_genes[[cluster]])))
    }
    
    wb <- createWorkbook()
    addWorksheet(wb, "top_genes")
    writeData(wb, "top_genes", top_genes_df)
    saveWorkbook(wb, file.path(output_dir, "top_genes.xlsx"), overwrite = TRUE)
  }
  
  # Create a summary file
  summary_text <- c(
    "Subset Analysis Results Summary",
    "============================",
    "",
    "Comparison Details:",
    paste("  Group 1:", analysis_results$metadata$group1_name, 
          sprintf("(n=%d cells)", analysis_results$metadata$n_cells_group1)),
    paste("  Group 2:", analysis_results$metadata$group2_name,
          sprintf("(n=%d cells)", analysis_results$metadata$n_cells_group2)),
    paste("  Mode:", analysis_results$metadata$mode),
    "",
    "Results:",
    paste("  Overall DE genes:", 
          ifelse(!is.null(analysis_results$overall_de), 
                 nrow(analysis_results$overall_de), 0)),
    paste("  Significant overall markers:", 
          ifelse(!is.null(analysis_results$overall_de_significant), 
                 nrow(analysis_results$overall_de_significant), 0)),
    paste("  Clusters analyzed:", length(analysis_results$cluster_specific)),
    paste("  Clusters with significant markers:",
          ifelse(!is.null(analysis_results$cluster_specific_significant),
                 length(analysis_results$cluster_specific_significant), 0)),
    paste("  Conserved markers:", 
          ifelse(!is.null(analysis_results$conserved_markers), 
                 nrow(analysis_results$conserved_markers), 0)),
    paste("  Significant conserved markers:", 
          ifelse(!is.null(analysis_results$conserved_markers_significant), 
                 nrow(analysis_results$conserved_markers_significant), 0)),
    "",
    "Files generated:",
    "- overall_DE.xlsx: Overall differential expression results (all genes)",
    paste0("    * Positive log2FC = higher in ", analysis_results$metadata$group1_name),
    paste0("    * Negative log2FC = higher in ", analysis_results$metadata$group2_name),
    "- overall_significant.xlsx: Significant overall markers (|avg_log2FC| > 0.5, p_val_adj < 0.01)",
    paste0("    * upregulated tab: Genes higher in ", analysis_results$metadata$group1_name),
    paste0("    * downregulated tab: Genes higher in ", analysis_results$metadata$group2_name),
    "- cluster_specific_DE.xlsx: Differential expression results for each cluster (all genes)",
    "- cluster_specific_significant.xlsx: Significant markers per cluster (|avg_log2FC| > 0.5, p_val_adj < 0.01)",
    "    * One tab per cluster and direction (e.g., cluster_0_up, cluster_0_down)",
    "    * Each tab sorted by fold change magnitude",
    "- conserved_markers.xlsx: Markers conserved across clusters",
    "- conserved_markers_significant.xlsx: Significant conserved markers (|avg_log2FC| > 0.5, p_val_adj < 0.01)",
    "    * upregulated tab: Genes with positive fold change (sorted high to low)",
    "    * downregulated tab: Genes with negative fold change (sorted low to high)",
    "- pseudobulk_expression.xlsx: Pseudobulk expression profiles",
    "- top_genes.xlsx: Top differential genes by cluster"
  )
  
  writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
  
  cat("Results saved to directory:", output_dir, "\n")
  cat("Generated files:\n")
  list.files(output_dir, pattern = "\\.xlsx$|\\.txt$") %>% 
    paste0("- ", .) %>%
    cat(sep = "\n")
}


# Saturation
## Function to perform sampling analysis on cellID data without reads
sample_cellID <- function(df, percentages = seq(10, 99, by = 10), replicates = 5) {
  # Get total number of rows and unique cellIDs in the full dataset
  total_rows <- nrow(df)
  total_unique_cellIDs <- length(unique(df$cellID))
  
  # Create a data frame to store results
  results <- data.frame(
    percentage = numeric(),
    replicate = numeric(),
    unique_cellIDs = numeric(),
    percent_coverage = numeric()
  )
  
  # Loop through each sampling percentage
  for (pct in percentages) {
    # Calculate sample size
    sample_size <- round(total_rows * (pct/100))
    
    # Perform multiple replicates for each percentage
    for (rep in 1:replicates) {
      # Sample rows with replacement
      sampled_indices <- sample(1:total_rows, size = sample_size, replace = TRUE)
      sampled_data <- df[sampled_indices, ]
      
      # Count unique cellIDs in the sample
      unique_cellIDs_sampled <- length(unique(sampled_data$cellID))
      
      # Calculate percentage of unique cellIDs covered
      percent_coverage <- (unique_cellIDs_sampled / total_unique_cellIDs) * 100
      
      # Add to results
      results <- rbind(results, data.frame(
        percentage = pct,
        replicate = rep,
        unique_cellIDs = unique_cellIDs_sampled,
        percent_coverage = percent_coverage
      ))
    }
  }
  
  return(results)
}
## Function to summarize and visualize results from sample_cellID
analyze_cellID_sampling <- function(results, title = "cellID Saturation Analysis") {
  # Calculate summary statistics
  summary_stats <- aggregate(
    percent_coverage ~ percentage, 
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  
  # Convert the summary statistics to a more readable format
  summary_df <- data.frame(
    percentage = summary_stats$percentage,
    mean_coverage = summary_stats$percent_coverage[, "mean"],
    sd_coverage = summary_stats$percent_coverage[, "sd"]
  )
  
  # Print the summary
  print(summary_df)
  
  # Plot the results if running in an environment with plotting capabilities
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    
    p <- ggplot(results, aes(x = percentage, y = percent_coverage)) +
      stat_summary(fun = mean, geom = "line", size = 1, color = "blue") +
      stat_summary(fun.data = function(x) {
        return(c(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x)))
      }, geom = "ribbon", alpha = 0.3) +
      geom_point(alpha = 0.2) +
      labs(
        title = title,
        x = "Sampling Percentage",
        y = "Percentage of Unique cellIDs Covered",
        caption = paste("Based on", length(unique(results$replicate)), "replicates per sampling percentage")
      ) +
      theme_minimal() +
      theme(aspect.ratio = 1) +  # Make plot area square
      scale_x_continuous(breaks = seq(0, 100, by = 10)) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))
    
    print(p)
  } else {
    message("Install ggplot2 package for visualization: install.packages('ggplot2')")
  }
  
  return(summary_df)
}
## Function to perform sampling analysis on UMI-cellID data with reads
sample_UMI_weighted <- function(df, percentages = seq(10, 99, by = 10), replicates = 5) {
  # Get total number of reads and unique cellID,UMI pairs in the full dataset
  total_reads <- sum(df$reads)
  total_unique_pairs <- nrow(df)
  total_unique_cellIDs <- length(unique(df$cellID))
  
  cat("Dataset info - Total reads:", total_reads, "Unique pairs:", total_unique_pairs, "\n")
  
  # Create a data frame to store results
  results <- data.frame(
    percentage = numeric(),
    replicate = numeric(),
    unique_pairs = numeric(),
    unique_cellIDs = numeric(),
    pair_coverage = numeric(),
    cellID_coverage = numeric()
  )
  
  # Pre-calculate weights once
  weights <- df$reads / total_reads
  
  # Loop through each sampling percentage
  for (pct in percentages) {
    cat("Processing", pct, "% sampling...\n")
    
    # Calculate sample size (number of reads to sample)
    sample_size <- round(total_reads * (pct/100))
    
    # If sample size is larger than total reads, just use all data
    if (sample_size >= total_reads) {
      for (rep in 1:replicates) {
        results <- rbind(results, data.frame(
          percentage = pct,
          replicate = rep,
          unique_pairs = total_unique_pairs,
          unique_cellIDs = total_unique_cellIDs,
          pair_coverage = 100,
          cellID_coverage = 100
        ))
      }
      next
    }
    
    # Perform multiple replicates for each percentage
    for (rep in 1:replicates) {
      start_time <- Sys.time()
      
      # More efficient sampling approach
      if (sample_size > nrow(df) * 10) {
        # For very large sample sizes, use a different approach
        # Sample based on cumulative probabilities
        cumprobs <- cumsum(weights)
        random_vals <- sort(runif(sample_size))
        sampled_indices <- findInterval(random_vals, cumprobs) + 1
      } else {
        # Use standard sampling for smaller sizes
        sampled_indices <- sample(1:nrow(df), size = sample_size, replace = TRUE, prob = weights)
      }
      
      sampled_data <- df[sampled_indices, ]
      
      # Count unique cellID,UMI pairs and cellIDs in the sample
      unique_pairs_sampled <- nrow(unique(sampled_data[, c("cellID", "UMI")]))
      unique_cellIDs_sampled <- length(unique(sampled_data$cellID))
      
      # Calculate coverage percentages
      pair_coverage <- (unique_pairs_sampled / total_unique_pairs) * 100
      cellID_coverage <- (unique_cellIDs_sampled / total_unique_cellIDs) * 100
      
      # Add to results
      results <- rbind(results, data.frame(
        percentage = pct,
        replicate = rep,
        unique_pairs = unique_pairs_sampled,
        unique_cellIDs = unique_cellIDs_sampled,
        pair_coverage = pair_coverage,
        cellID_coverage = cellID_coverage
      ))
      
      # elapsed <- Sys.time() - start_time
      # if (elapsed > 5) {  # Warn if taking more than 5 seconds per replicate
      #   cat("  Replicate", rep, "took", round(elapsed, 1), "seconds\n")
      # }
    }
  }
  
  return(results)
}
##Function to summarize and visualize results with both pair and cellID coverage
analyze_UMI_sampling <- function(results, title = "Sampling Analysis with Read Weights") {
  # Calculate summary statistics for both pair and cellID coverage
  pair_summary <- aggregate(
    pair_coverage ~ percentage, 
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  
  cellID_summary <- aggregate(
    cellID_coverage ~ percentage, 
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  
  # Convert the summary statistics to a more readable format
  summary_df <- data.frame(
    percentage = pair_summary$percentage,
    mean_pair_coverage = pair_summary$pair_coverage[, "mean"],
    sd_pair_coverage = pair_summary$pair_coverage[, "sd"],
    mean_cellID_coverage = cellID_summary$cellID_coverage[, "mean"],
    sd_cellID_coverage = cellID_summary$cellID_coverage[, "sd"]
  )
  
  # Print the summary
  print(summary_df)
  
  # Plot the results if running in an environment with plotting capabilities
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    
    # Reshape data for better plotting
    plot_data <- rbind(
      data.frame(results, type = "cellID Coverage", coverage = results$cellID_coverage),
      data.frame(results, type = "cellID,UMI Pair Coverage", coverage = results$pair_coverage)
    )
    
    p <- ggplot(plot_data, aes(x = percentage, y = coverage, color = type, fill = type)) +
      stat_summary(fun = mean, geom = "line", size = 1, aes(group = type)) +
      stat_summary(fun.data = function(x) {
        return(c(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x)))
      }, geom = "ribbon", alpha = 0.2, aes(group = type)) +
      labs(
        title = title,
        x = "Percentage of Total Reads Sampled",
        y = "Coverage (%)",
        caption = paste("Based on", length(unique(results$replicate)), "replicates per sampling percentage")
      ) +
      theme_minimal() +
      theme(
        legend.title = element_blank(),
        aspect.ratio = 0.8
      ) +
      theme(aspect.ratio = 1) +  # Make plot area square
      scale_x_continuous(breaks = seq(0, 100, by = 10)) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))
    
    print(p)
  } else {
    message("Install ggplot2 package for visualization: install.packages('ggplot2')")
  }
  
  return(summary_df)
}
## Function to fit models to results of sample_UMI
fit_models <- function(results, target_coverage = 100, coverage_column = NULL) {
  # Check the structure of the results dataframe to determine which column to use
  if (is.null(coverage_column)) {
    if ("percent_coverage" %in% colnames(results)) {
      coverage_column <- "percent_coverage"
      message("Using 'percent_coverage' column")
    } else if ("pair_coverage" %in% colnames(results)) {
      coverage_column <- "pair_coverage"
      message("Using 'pair_coverage' column")
    } else if ("cellID_coverage" %in% colnames(results)) {
      coverage_column <- "cellID_coverage"
      message("Using 'cellID_coverage' column")
    } else {
      stop("No coverage column found in results. Please specify 'coverage_column'")
    }
  } else {
    if (!coverage_column %in% colnames(results)) {
      stop(paste("Specified column", coverage_column, "not found in results dataframe"))
    }
  }
  
  # Create a formula for aggregation
  agg_formula <- as.formula(paste(coverage_column, "~ percentage"))
  
  # Aggregate to get mean coverage by percentage
  agg_data <- aggregate(agg_formula, data = results, FUN = mean)
  
  # Rename column for consistent processing
  names(agg_data)[names(agg_data) == coverage_column] <- "coverage"
  
  # Remove any rows with invalid coverage values
  agg_data <- agg_data[is.finite(agg_data$coverage) & agg_data$coverage >= 0, ]
  
  if (nrow(agg_data) == 0) {
    stop("No valid coverage data found after cleaning")
  }
  
  # Fit several models
  models <- list()
  
  # 1. Michaelis-Menten model (y = (Vmax * x) / (Km + x))
  mm_model <- try({
    nls(coverage ~ (Vmax * percentage) / (Km + percentage),
        data = agg_data,
        start = list(Vmax = max(agg_data$coverage) * 1.1, 
                     Km = median(agg_data$percentage)),
        control = nls.control(maxiter = 200))
  }, silent = TRUE)
  
  if (!inherits(mm_model, "try-error")) {
    models$michaelis_menten <- mm_model
  }
  
  # 2. Log model (y = a + b * log(x)) - only if percentage > 0
  if (min(agg_data$percentage) > 0) {
    log_model <- try({
      nls(coverage ~ a + b * log(percentage),
          data = agg_data,
          start = list(a = min(agg_data$coverage), 
                       b = (max(agg_data$coverage) - min(agg_data$coverage)) / 
                         log(max(agg_data$percentage))),
          control = nls.control(maxiter = 200))
    }, silent = TRUE)
    
    if (!inherits(log_model, "try-error")) {
      models$logarithmic <- log_model
    }
  }
  
  # 3. Power model (y = a * x^b)
  if (min(agg_data$percentage) > 0 && min(agg_data$coverage) > 0) {
    power_model <- try({
      nls(coverage ~ a * percentage^b,
          data = agg_data,
          start = list(a = 1, b = 0.5),
          control = nls.control(maxiter = 200))
    }, silent = TRUE)
    
    if (!inherits(power_model, "try-error")) {
      models$power <- power_model
    }
  }
  
  # If all models failed, try a simple loess smoother
  if (length(models) == 0) {
    message("All nonlinear models failed to converge. Using loess smoothing instead.")
    loess_model <- loess(coverage ~ percentage, data = agg_data, span = 0.75)
    models$loess <- loess_model
  }
  
  # Create prediction grid for plotting
  pred_grid <- data.frame(percentage = seq(max(1, min(agg_data$percentage)), 
                                           max(agg_data$percentage) * 1.2, 
                                           length.out = 100))
  
  # Make predictions for each model
  for (model_name in names(models)) {
    predicted_values <- try({
      if (model_name == "loess") {
        predict(models[[model_name]], newdata = pred_grid)
      } else {
        predict(models[[model_name]], newdata = pred_grid)
      }
    }, silent = TRUE)
    
    if (!inherits(predicted_values, "try-error")) {
      # Clean predictions - remove negative values and extreme outliers
      predicted_values[predicted_values < 0] <- 0
      predicted_values[predicted_values > 200] <- NA  # Remove unrealistic values
      pred_grid[[model_name]] <- predicted_values
    }
  }
  
  # Create a plot with all models
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    library(reshape2)
    
    # Only include columns that actually have predictions
    pred_cols <- intersect(names(models), colnames(pred_grid))
    if (length(pred_cols) > 0) {
      pred_subset <- pred_grid[, c("percentage", pred_cols)]
      
      # Melt the prediction data for plotting
      pred_long <- melt(pred_subset, id.vars = "percentage", 
                        variable.name = "model", value.name = "predicted_coverage")
      
      # Remove rows with NA predictions
      pred_long <- pred_long[!is.na(pred_long$predicted_coverage), ]
      
      if (nrow(pred_long) > 0) {
        # Get title based on coverage column used
        plot_title <- paste("Model Fitting for", gsub("_", " ", coverage_column))
        y_label <- paste("Percentage of", gsub("_", " ", coverage_column))
        
        # Create the plot
        p <- ggplot() +
          geom_point(data = agg_data, aes(x = percentage, y = coverage), 
                     alpha = 0.6, size = 3) +
          geom_line(data = pred_long, aes(x = percentage, y = predicted_coverage, 
                                          color = model), size = 1) +
          geom_hline(yintercept = target_coverage, linetype = "dashed", color = "red") +
          labs(
            title = plot_title,
            x = "Sampling Percentage",
            y = y_label,
            color = "Model"
          ) +
          theme_minimal() +
          theme(aspect.ratio = 1) +
          scale_x_continuous(breaks = seq(0, max(pred_grid$percentage), by = 20)) +
          scale_y_continuous(limits = c(0, min(105, max(c(agg_data$coverage, 
                                                          pred_long$predicted_coverage), 
                                                        na.rm = TRUE) * 1.1)), 
                             breaks = seq(0, 100, by = 10))
        
        print(p)
      }
    }
  }
  
  # Estimate required percentage for target coverage
  result_table <- data.frame(
    model = character(),
    required_percentage = numeric(),
    max_achievable = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (model_name in names(models)) {
    if (model_name %in% colnames(pred_grid)) {
      predicted <- pred_grid[[model_name]]
      predicted <- predicted[is.finite(predicted)]
      
      if (length(predicted) > 0) {
        # Find max achievable coverage
        max_achievable <- max(predicted, na.rm = TRUE)
        
        # Calculate required percentage
        if (max_achievable >= target_coverage) {
          # Find the first percentage where we achieve target coverage
          target_indices <- which(predicted >= target_coverage)
          if (length(target_indices) > 0) {
            required_pct <- pred_grid$percentage[min(target_indices)]
          } else {
            required_pct <- NA
          }
        } else {
          required_pct <- NA
        }
        
        # Add to results
        result_table <- rbind(result_table, data.frame(
          model = model_name,
          required_percentage = required_pct,
          max_achievable = max_achievable
        ))
      }
    }
  }
  
  # Print results
  if (nrow(result_table) > 0) {
    print(result_table)
  } else {
    message("No valid model predictions available")
  }
  
  return(list(
    models = models, 
    predictions = pred_grid, 
    required = result_table, 
    coverage_column = coverage_column,
    aggregated_data = agg_data
  ))
}


### Prognostic functions: let Claude tell you what to do!
analyze_sequencing_value <- function(model_results, target_coverage = 90) {
  
  # Get coverage data from successful model predictions
  pred_data <- model_results$predictions
  
  # Find all model prediction columns (exclude 'percentage')
  model_columns <- setdiff(colnames(pred_data), "percentage")
  valid_maxes <- c()
  
  # Extract maximum coverage from each valid model
  for (col in model_columns) {
    if (col %in% colnames(pred_data)) {
      model_data <- pred_data[[col]]
      # Remove infinite, NA, and unrealistic values
      model_data <- model_data[is.finite(model_data) & model_data >= 0 & model_data <= 200]
      
      if (length(model_data) > 0) {
        model_max <- max(model_data, na.rm = TRUE)
        valid_maxes <- c(valid_maxes, model_max)
      }
    }
  }
  
  # Use the best estimate from valid models
  if (length(valid_maxes) > 0) {
    current_max_coverage <- max(valid_maxes)
  } else {
    # Fallback: try to get coverage from aggregated data
    if ("aggregated_data" %in% names(model_results) && 
        nrow(model_results$aggregated_data) > 0) {
      current_max_coverage <- max(model_results$aggregated_data$coverage, na.rm = TRUE)
      message("Using maximum observed coverage from data since models failed")
    } else {
      current_max_coverage <- NA
    }
  }
  
  cat(sprintf("Overall max coverage from valid models: %.1f%%\n", current_max_coverage))
  
  # Calculate and display saturation levels
  if (!is.na(current_max_coverage)) {
    sequencing_saturation <- min(100, current_max_coverage)  # Current sequencing relative to max possible
    cat(sprintf("📊 SATURATION ANALYSIS:\n"))
    cat(sprintf("   Sequencing saturation: %.1f%% (how much of library complexity you've captured)\n", sequencing_saturation))
    cat(sprintf("   Library saturation: %.1f%% (max coverage this library can achieve)\n", current_max_coverage))
    cat("\n")
  }
  
  # Check if target is achievable
  if (is.na(current_max_coverage) || current_max_coverage < target_coverage) {
    cat("❌ DECISION A (More sequencing): NOT RECOMMENDED\n")
    if (is.na(current_max_coverage)) {
      cat("   Unable to determine library coverage (model fitting failed)\n")
    } else {
      cat(sprintf("   Your current library can only achieve %.1f%% coverage\n", current_max_coverage))
      cat(sprintf("   Even infinite sequencing won't reach your %.0f%% target\n", target_coverage))
    }
    cat("   -> This library is molecularly limited\n\n")
    
    return(list(
      more_sequencing_recommended = FALSE,
      reason = ifelse(is.na(current_max_coverage), "Model fitting failed", "Library complexity limitation"),
      max_achievable = current_max_coverage,
      already_sufficient = FALSE
    ))
  }
  
  # If we can achieve the target, find required sampling percentage
  cat("✅ Library can achieve target coverage!\n")
  
  # Find required sampling percentage from results table
  required_results <- model_results$required
  valid_requirements <- required_results[is.finite(required_results$required_percentage) & 
                                           !is.na(required_results$required_percentage), ]
  
  if (nrow(valid_requirements) > 0) {
    # Take the most conservative (highest) estimate
    required_pct <- max(valid_requirements$required_percentage, na.rm = TRUE)
    
    # Calculate percentage increase needed (not fold increase)
    percent_increase <- required_pct - 100
    
    if (required_pct <= 100) {
      cat("✅ DECISION A (More sequencing): RECOMMENDED\n")
      cat(sprintf("   You need ~%.0f%% sampling to reach %.0f%% coverage\n", 
                  required_pct, target_coverage))
      if (required_pct < 100) {
        cat(sprintf("   You're already sequencing enough (could reduce by %.0f%%)\n", 
                    100 - required_pct))
      } else {
        cat("   Your current sequencing depth is optimal\n")
      }
    } else if (required_pct <= 120) {
      cat("⚠️  DECISION A (More sequencing): MARGINAL VALUE\n")
      cat(sprintf("   You need %.0f%% sampling to reach %.0f%% coverage\n", 
                  required_pct, target_coverage))
      cat(sprintf("   Increase sequencing by %.0f%% more reads\n", percent_increase))
      cat("   This is close to your library's complexity limit\n")
    } else {
      cat("❌ DECISION A (More sequencing): NOT RECOMMENDED\n")
      cat(sprintf("   You would need %.0f%% sampling to reach %.0f%% coverage\n", 
                  required_pct, target_coverage))
      cat(sprintf("   This would require %.0f%% more reads\n", percent_increase))
      cat("   This indicates you're near the library complexity limit\n")
    }
    
    # Check if already sufficient
    already_sufficient <- required_pct <= 100
    
  } else {
    required_pct <- NA
    percent_increase <- NA
    already_sufficient <- FALSE
    cat("⚠️  DECISION A (More sequencing): UNCERTAIN\n")
    cat("   Models couldn't determine required sampling percentage\n")
    cat("   But your library appears capable of reaching the target\n")
  }
  
  return(list(
    more_sequencing_recommended = !is.na(required_pct) && required_pct <= 120,
    required_fold_increase = ifelse(is.na(required_pct), NA, required_pct/100),
    percent_increase_needed = percent_increase,
    max_achievable = current_max_coverage,
    already_sufficient = already_sufficient
  ))
}

analyze_library_scaling <- function(reads_GEX, model_results, target_coverage = 90) {
  
  # Current library stats
  current_cells <- length(unique(reads_GEX$cellID))
  current_pairs <- nrow(reads_GEX)
  current_reads <- sum(reads_GEX$reads)
  
  # Efficiency metrics
  pairs_per_cell <- current_pairs / current_cells
  reads_per_cell <- current_reads / current_cells
  
  # Get max coverage safely from model predictions
  pred_data <- model_results$predictions
  model_columns <- setdiff(colnames(pred_data), "percentage")
  valid_maxes <- c()
  
  for (col in model_columns) {
    if (col %in% colnames(pred_data)) {
      model_data <- pred_data[[col]]
      # Remove infinite, NA, and negative values
      model_data <- model_data[is.finite(model_data) & model_data >= 0]
      
      if (length(model_data) > 0) {
        model_max <- max(model_data, na.rm = TRUE)
        valid_maxes <- c(valid_maxes, model_max)
      }
    }
  }
  
  # Take the best estimate from valid models
  if (length(valid_maxes) > 0) {
    current_max <- max(valid_maxes)
  } else {
    current_max <- NA
  }
  
  cat("📊 CURRENT LIBRARY EFFICIENCY:\n")
  cat(sprintf("   Cells captured: %d\n", current_cells))
  cat(sprintf("   Unique pairs per cell: %.1f\n", pairs_per_cell))
  cat(sprintf("   Reads per cell: %.1f\n", reads_per_cell))
  if (!is.na(current_max)) {
    cat(sprintf("   Current max coverage: %.1f%%\n", current_max))
  } else {
    cat("   Current max coverage: Unable to determine\n")
  }
  
  # Estimate scaling needed
  if (is.na(current_max)) {
    cat("\n⚠️  DECISION B (Bigger library): UNCERTAIN\n")
    cat("   Unable to determine current coverage efficiency\n")
    recommend_bigger <- TRUE  # Conservative recommendation
    scaling_factor <- 2  # Default guess
  } else if (current_max >= target_coverage) {
    cat("\n✅ DECISION B (Bigger library): NOT NEEDED\n")
    cat("   Your current library size is sufficient\n")
    recommend_bigger <- FALSE
    scaling_factor <- 1
  } else {
    # Rough estimate: if we're capturing X% of cells efficiently,
    # we might need ~1/X times more input to capture target%
    scaling_factor <- target_coverage / current_max
    
    cat(sprintf("\n🔄 DECISION B (Bigger library): RECOMMENDED\n"))
    cat(sprintf("   Estimated scaling needed: %.1fx more input\n", scaling_factor))
    cat(sprintf("   Try ~%.0f cells instead of %d\n", 
                current_cells * scaling_factor, current_cells))
    recommend_bigger <- TRUE
  }
  
  return(list(
    bigger_library_recommended = recommend_bigger,
    current_efficiency = list(
      cells = current_cells,
      pairs_per_cell = pairs_per_cell,
      max_coverage = current_max
    ),
    scaling_estimate = scaling_factor
  ))
}

make_sequencing_decision <- function(reads_GEX, target_coverage = 90, 
                                     percentages = c(seq(10, 80, by = 10), 
                                                     seq(85, 99, by = 2.5)),
                                     coverage_column = NULL) {
  
  tryCatch({
    results <- sample_UMI_weighted(reads_GEX, percentages = percentages)
    
    # Check if results are valid
    if (is.null(results) || nrow(results) == 0) {
      stop("No results from sampling function")
    }
    
    model_results <- fit_models(results, target_coverage = target_coverage, coverage_column = coverage_column)
    
  }, error = function(e) {
    cat("Error in analysis pipeline:\n")
    cat(paste("  ", e$message, "\n"))
    return(NULL)
  })
  
  # Analyze both scenarios
  seq_analysis <- analyze_sequencing_value(model_results, target_coverage)
  lib_analysis <- analyze_library_scaling(reads_GEX, model_results, target_coverage)
  
  # Print summary with fixed string concatenation
  separator <- paste(rep("=", 50), collapse = "")
  cat(paste(separator, "\n", sep = ""))
  
  # Decision tree
  seq_recommended <- seq_analysis$more_sequencing_recommended
  lib_recommended <- lib_analysis$bigger_library_recommended
  already_sufficient <- seq_analysis$already_sufficient
  
  if (!is.na(already_sufficient) && already_sufficient) {
    cat("✅ RECOMMENDATION: Current sequencing depth is sufficient\n")
    cat(sprintf("   You have already achieved/exceeded the %.0f%% coverage target\n", target_coverage))
    if (!is.na(seq_analysis$required_fold_increase) && seq_analysis$required_fold_increase < 0.8) {
      cat("   You may be sequencing to excess - consider reducing depth for cost savings\n")
    }
    
  } else if (!is.na(seq_recommended) && seq_recommended && !lib_recommended) {
    cat("📈 RECOMMENDATION: Sequence deeper\n")
    if (!is.na(seq_analysis$percent_increase_needed)) {
      cat(sprintf("   Increase sequencing by %.0f%% more reads to reach %.0f%% coverage\n",
                  seq_analysis$percent_increase_needed, target_coverage))
    } else if (!is.na(seq_analysis$required_fold_increase)) {
      # Fallback to fold increase if percent increase not available
      percent_inc <- (seq_analysis$required_fold_increase - 1) * 100
      cat(sprintf("   Increase sequencing by %.0f%% more reads to reach %.0f%% coverage\n",
                  percent_inc, target_coverage))
    }
    
  } else if ((is.na(seq_recommended) || !seq_recommended) && lib_recommended) {
    cat("🔬 RECOMMENDATION: Make bigger library\n")
    if (!is.na(seq_analysis$max_achievable)) {
      cat(sprintf("   Your current library is complexity-limited at %.1f%%\n",
                  seq_analysis$max_achievable))
    }
    scaling_percent <- (lib_analysis$scaling_estimate - 1) * 100
    cat(sprintf("   Try %.0f%% more input material\n", scaling_percent))
    
  } else if (!is.na(seq_recommended) && seq_recommended && lib_recommended) {
    cat("⚖️  RECOMMENDATION: Both options viable\n")
    cat("   Choice depends on cost and timeline:\n")
    if (!is.na(seq_analysis$percent_increase_needed)) {
      cat(sprintf("   - More sequencing: %.0f%% more reads\n", 
                  seq_analysis$percent_increase_needed))
    }
    scaling_percent <- (lib_analysis$scaling_estimate - 1) * 100
    cat(sprintf("   - Bigger library: %.0f%% more input material\n", scaling_percent))
    
  } else {
    cat("✅ RECOMMENDATION: Current setup is sufficient\n")
    if (!is.na(seq_analysis$max_achievable)) {
      cat(sprintf("   You can already achieve %.1f%% coverage\n",
                  seq_analysis$max_achievable))
    }
  }
  
  return(list(
    sequencing = seq_analysis,
    library = lib_analysis,
    target_coverage = target_coverage,
    model_results = model_results
  ))
}
