# functions.R

# Loading
## standard processing of counts matrices with Seurat
SeuratPipeline <- function(file_location, sample_name, output_dir = "/projects/b1042/GoyalLab/egrody/", plots = FALSE, rds = FALSE) {
  #if you would like to output plots from this pipeline, pass entries for the plots and output_dir variables
  
  # Reading in feature matrices
  data <- Read10X(data.dir = file_location)
  # Initialize the Seurat object with the raw (non-normalized data)
  sample <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  # QC
  sample <- Add_Mito_Ribo(sample, species = "macaque")
  if (plots) {
    vplot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
    ggsave(vplot, file = paste0(output_dir, sample_name, "_preQC_violinplot.svg"))
  }
  sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA < 20000 & percent_mito < 5)
  if (plots) {
    vplot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
    ggsave(vplot, file = paste0(output_dir, sample_name, "_postQC_violinplot.svg"))
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
  sample <- FindClusters(sample, resolution = 0.5, verbose = FALSE)
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
## processing SALVEseq count matrices with Seurat
SeuratSALVEseq <- function(file_location, sample_name, output_dir = "/projects/b1042/GoyalLab/egrody/", plots = FALSE, rds = FALSE) {
  #if you would like to output plots from this pipeline, pass entries for the plots and output_dir variables
  
  # Reading in feature matrices
  data <- Read10X(data.dir = file_location)
  # Initialize the Seurat object with the raw (non-normalized data)
  sample <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3) #no minimum features 
  # QC
  sample <- Add_Mito_Ribo(sample, species = "macaque")
  if (plots) {
    vplot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
    ggsave(vplot, file = paste0(output_dir, sample_name, "_noQC_violinplot.svg"))
  }
  # Normalizing
  sample <- NormalizeData(sample, verbose = FALSE)
  # Return the variable
  return(sample)
}

# Output Dataframes
## isolate specific gene(s) expression from a Seurat object
targetExpressionDF <- function(seurat_obj, genes, count_type = "normalized") {
  if (!is.vector(genes)) {
    genes <- c(genes)
  }
  
  if (count_type == "raw") {
    expression <- GetAssayData(object = seurat_obj, assay = "RNA", layer = "counts")[genes,, drop=FALSE]
  } else {
    expression <- GetAssayData(object = seurat_obj, assay = "RNA", layer = "data")[genes,, drop=FALSE]
  }
  
  if (length(genes) == 1) {
    expression_df <- data.frame(
      cellID = colnames(expression),
      setNames(list(as.numeric(expression)), genes),
      stringsAsFactors = FALSE
    )
  } else {
    expression_df <- data.frame(
      cellID = colnames(expression),
      t(expression),
      stringsAsFactors = FALSE
    )
  }
  
  # Check if UMAP reduction exists
  if ("umap" %in% names(seurat_obj@reductions)) {
    umap_coords <- Embeddings(object = seurat_obj, reduction = "umap")
    cluster <- Idents(object = seurat_obj)
    
    SingleCell <- expression_df %>%
      mutate(
        UMAP1 = umap_coords[, 1],
        UMAP2 = umap_coords[, 2],
        cluster = cluster
      )
  } else {
    warning("UMAP not present, generating dataframe without UMAP coordinates.")
    SingleCell <- expression_df
  }
  
  return(SingleCell)
}
## load scPathoQuant data
scPathoCountCellID <- function(seurat_obj, viral_gene, count_type = "normalized") {
  # Get expression data
  if (count_type == "raw") {
    expression <- GetAssayData(object = seurat_obj, assay = "RNA", slot = "counts")[viral_gene,]
  } else {
    expression <- GetAssayData(object = seurat_obj, assay = "RNA", slot = "data")[viral_gene,]
  }
  # Get UMAP coordinates
  umap_coords <- Embeddings(object = seurat_obj, reduction = "umap")
  
  SingleCell <- data.frame(
    cellID = names(expression),
    Count = unname(expression),
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    stringsAsFactors = FALSE
  )
  
  return(SingleCell)
}
## output viral expression .csvs
save_scPQ_expression <- function(seurat_objects, gene_names, output_dir = ".") {
  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (obj_name in names(seurat_objects)) {
    seurat_obj <- seurat_objects[[obj_name]]
    
    # Check which genes are present in the Seurat object
    genes_present <- gene_names[gene_names %in% rownames(seurat_obj)]
    genes_missing <- setdiff(gene_names, genes_present)
    
    if (length(genes_missing) > 0) {
      warning(paste("In object", obj_name, "the following genes are missing:", 
                    paste(genes_missing, collapse = ", "), "\n"))
    }
    
    for (gene in genes_present) {
      tryCatch({
        # Run scPathoCountCellID function
        result <- scPathoCountCellID(seurat_obj, gene)
        
        file_path <- paste0(output_dir, obj_name, "_scPQexpression_", gene, ".csv")
        write.csv(result, file_path)
        
        cat("Saved:", file_path, "\n")
      }, error = function(e) {
        warning(paste("Error processing gene", gene, "in object", 
                      obj_name, ":", e$message, "/n"))
      })
    }
  }
}

# Joint Analyses
## joining and printing numbers

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

# Plotting
## plot a UMAP
plotUMAP <- function(data, colorby, title, output_dir, saveas, comparison = FALSE) {
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
    umap <- umap + scale_color_gradient(low = "lightgrey", high = "darkblue")
  }
  ggsave(umap, file = paste0(output_dir, saveas))
  
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
barcodeRankPlot <- function(rawDataFolder, jointFullJoin, plotTitle = "Barcode Rank Plot", output_dir = "/projects/b1042/GoyalLab/egrody/", saveas = "barcodeRankPlot.svg") {
  # add error for wrong format of jointFullJoin columns
  if (!"log1pSALVE" %in% colnames(jointFullJoin)) {
    stop("Error: jointFullJoin must contain a 'log1pSALVE' column")
  }
  
  if (!"cellID" %in% colnames(jointFullJoin)) {
    stop("Error: jointFullJoin must contain a 'cellID' column")
  }
  # Create a data frame with barcode and UMI count
  rawdata <- Read10X(rawDataFolder)
  umi_counts <- colSums(rawdata)
  plot_data <- data.frame(barcode = names(umi_counts), 
                          umi_count = umi_counts)
  plot_data <- plot_data %>% 
    arrange(desc(umi_count)) %>%
    mutate(barcode = substr(barcode, 1, nchar(barcode) - 2))
  
  # highlight the SALVE+ cellIDs
  SALVE_cellids <- jointFullJoin %>% filter(log1pSALVE > 0) %>% select(cellID)
  cat(sum(SALVE_cellids[[1]] %in% plot_data$barcode), "out of", nrow(SALVE_cellids), 
      "SALVEseq cellIDs are found in raw cellID list\n")
  plot_data$highlight <- plot_data$barcode %in% SALVE_cellids[[1]]
  
  #count SALVE+ cellIDs with transcriptomes
  SALVE_in_rawdata <- plot_data %>% 
    filter(barcode %in% SALVE_cellids[[1]])
  
  num_full_transcriptomes <- sum(SALVE_in_rawdata$umi_count >= 100)
  cat(num_full_transcriptomes, "out of", nrow(SALVE_cellids), 
      "SALVEseq cellIDs have full transcriptomes (>= 100 UMI)")
  
  set.seed(123)
  plot_data <- plot_data %>%
    arrange(desc(umi_count)) %>%
    mutate(
      rank = row_number(),
      umi_count_adj = ifelse(umi_count <= 0, 0.1, umi_count)  # Replace 0 or negative with 0.1
    ) %>%
    slice_sample(n = 20000) # Update if necessary
  
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
  ggsave(paste0(output_dir, saveas), umi_rank_plot, width = 10, height = 8, device = "svg")
  
  return(umi_rank_plot)
}

# DEG Analysis
## perform comprehensive subset analysis
analyze_cell_subset <- function(seurat_obj, subset_cells, min_pct = 0.1) {
  require(Seurat)
  require(dplyr)
  
  # Input validation
  if (!is(seurat_obj, "Seurat")) {
    stop("First argument must be a Seurat object")
  }
  
  if (!is.vector(subset_cells)) {
    stop("subset_cells must be a vector of cell IDs")
  }
  
  # Convert cell names to character vector if they aren't already
  subset_cells <- as.character(subset_cells)
  
  # Check if the subset cells exist in the Seurat object
  valid_cells <- subset_cells %in% colnames(seurat_obj)
  if (!any(valid_cells)) {
    stop("None of the provided cell IDs were found in the Seurat object")
  }
  
  if (!all(valid_cells)) {
    warning(sprintf("%d cells from the subset were not found in the Seurat object",
                    sum(!valid_cells)))
    subset_cells <- subset_cells[valid_cells]
  }
  
  # Print diagnostic information
  cat(sprintf("Found %d cells in the subset\n", length(subset_cells)))
  
  # Add metadata column for subset membership
  seurat_obj$in_subset <- colnames(seurat_obj) %in% subset_cells
  
  # Get unique clusters containing subset cells
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    stop("No 'seurat_clusters' column found in metadata. Please run FindClusters() first")
  }
  
  subset_clusters <- unique(seurat_obj$seurat_clusters[colnames(seurat_obj) %in% subset_cells])
  
  # Print cluster information
  cat(sprintf("Subset cells are found in %d clusters: %s\n",
              length(subset_clusters),
              paste(subset_clusters, collapse = ", ")))
  
  # Initialize results list
  results <- list()
  
  # 1. Cluster-specific analysis
  cluster_results <- list()
  for(cluster in subset_clusters) {
    # Create logical vectors for the comparison
    cluster_cells <- seurat_obj$seurat_clusters == cluster
    subset_in_cluster <- colnames(seurat_obj)[cluster_cells] %in% subset_cells
    
    # Check if we have enough cells in both groups for comparison
    n_subset <- sum(subset_in_cluster)
    n_other <- sum(cluster_cells) - n_subset
    
    if (n_subset < 3 || n_other < 3) {
      warning(sprintf("Skipping cluster %s: insufficient cells (subset: %d, other: %d)",
                      cluster, n_subset, n_other))
      next
    }
    
    # Create a temporary Seurat object for this cluster
    cluster_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cluster_cells])
    Idents(cluster_obj) <- factor(cluster_obj$in_subset)
    
    # Perform DE analysis within this cluster
    tryCatch({
      de_results <- FindMarkers(cluster_obj,
                                ident.1 = TRUE,
                                ident.2 = FALSE,
                                min.pct = min_pct,
                                test.use = "wilcox")
      
      cluster_results[[paste0("cluster_", cluster)]] <- de_results
      
      cat(sprintf("Successfully analyzed cluster %s (subset: %d, other: %d)\n",
                  cluster, n_subset, n_other))
      
    }, error = function(e) {
      warning(sprintf("Error in DE analysis for cluster %s: %s", cluster, e$message))
    })
  }
  results$cluster_specific <- cluster_results
  
  # 2. Calculate average expression profiles
  results$avg_expression <- tryCatch({
    AverageExpression(seurat_obj,
                      group.by = "in_subset",
                      assays = "RNA")$RNA
  }, error = function(e) {
    warning("Error calculating average expression: ", e$message)
    return(NULL)
  })
  
  # 3. Modified conserved markers analysis
  if (length(cluster_results) >= 2) {  # Only try if we have at least 2 clusters
    # Set identities for conserved marker analysis
    Idents(seurat_obj) <- factor(seurat_obj$in_subset)
    
    results$conserved_markers <- tryCatch({
      FindMarkers(seurat_obj,
                  ident.1 = TRUE,
                  ident.2 = FALSE,
                  min.pct = min_pct,
                  test.use = "wilcox",
                  subset.ident = subset_clusters)
    }, error = function(e) {
      warning("Error finding conserved markers: ", e$message)
      return(NULL)
    })
  } else {
    warning("Not enough clusters with sufficient cells for conserved marker analysis")
    results$conserved_markers <- NULL
  }
  
  # 4. Get top genes from differential expression results
  get_top_genes <- function(de_results, n = 40) {
    if (is.null(de_results) || nrow(de_results) == 0) return(character(0))
    top_genes <- rownames(de_results[de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 0.5, ])
    return(head(top_genes, n))
  }
  
  top_genes_by_cluster <- lapply(cluster_results, get_top_genes)
  results$top_genes <- top_genes_by_cluster
  
  # Summary of analysis
  cat("\nAnalysis Summary:\n")
  cat(sprintf("- Analyzed %d clusters\n", length(cluster_results)))
  cat(sprintf("- Found average expression profiles: %s\n", !is.null(results$avg_expression)))
  cat(sprintf("- Found conserved markers: %s\n", !is.null(results$conserved_markers)))
  
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
          height = h,
          dpi = dpi
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
  
  # 4. Save top genes list
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
    paste("Number of clusters analyzed:", length(analysis_results$cluster_specific)),
    paste("Number of conserved markers:", 
          ifelse(!is.null(analysis_results$conserved_markers), 
                 nrow(analysis_results$conserved_markers), 0)),
    "\nFiles generated:",
    "- cluster_specific_DE.xlsx: Differential expression results for each cluster",
    "- conserved_markers.xlsx: Markers conserved across clusters",
    "- average_expression.xlsx: Average expression profiles",
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
sample_cellID <- function(df, percentages = seq(10, 100, by = 10), replicates = 5) {
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
sample_UMI_weighted <- function(df, percentages = seq(10, 100, by = 10), replicates = 5) {
  # Get total number of reads and unique cellID,UMI pairs in the full dataset
  total_reads <- sum(df$reads)
  total_unique_pairs <- nrow(df)
  total_unique_cellIDs <- length(unique(df$cellID))
  
  # Create a data frame to store results
  results <- data.frame(
    percentage = numeric(),
    replicate = numeric(),
    unique_pairs = numeric(),
    unique_cellIDs = numeric(),
    pair_coverage = numeric(),
    cellID_coverage = numeric()
  )
  
  # Loop through each sampling percentage
  for (pct in percentages) {
    # Calculate sample size (number of reads to sample)
    sample_size <- round(total_reads * (pct/100))
    
    # Perform multiple replicates for each percentage
    for (rep in 1:replicates) {
      # Create probability weights based on read counts
      weights <- df$reads / total_reads
      
      # Sample rows with replacement according to read count weights
      sampled_indices <- sample(1:nrow(df), size = sample_size, replace = TRUE, prob = weights)
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
    } else if ("cellID_coverage" %in% colnames(results)) {
      coverage_column <- "cellID_coverage"
      message("Using 'cellID_coverage' column")
    } else if ("pair_coverage" %in% colnames(results)) {
      coverage_column <- "pair_coverage"
      message("Using 'pair_coverage' column")
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
  
  # Fit several models
  models <- list()
  
  # 1. Asymptotic model (y = a * (1 - exp(-b * x)))
  asymp_model <- try(nls(coverage ~ a * (1 - exp(-b * percentage)), 
                         data = agg_data,
                         start = list(a = 100, b = 0.05),
                         control = nls.control(maxiter = 200)),
                     silent = TRUE)
  
  if (!inherits(asymp_model, "try-error")) {
    models$asymptotic <- asymp_model
  }
  
  # 2. Michaelis-Menten model (y = (Vmax * x) / (Km + x))
  mm_model <- try(nls(coverage ~ (Vmax * percentage) / (Km + percentage),
                      data = agg_data,
                      start = list(Vmax = 100, Km = 20),
                      control = nls.control(maxiter = 200)),
                  silent = TRUE)
  
  if (!inherits(mm_model, "try-error")) {
    models$michaelis_menten <- mm_model
  }
  
  # 3. Log model (y = a + b * log(x))
  log_model <- try(nls(coverage ~ a + b * log(percentage),
                       data = agg_data,
                       start = list(a = 0, b = 20),
                       control = nls.control(maxiter = 200)),
                   silent = TRUE)
  
  if (!inherits(log_model, "try-error")) {
    models$logarithmic <- log_model
  }
  
  # 4. Logistic model (y = L / (1 + exp(-k * (x - x0))))
  logistic_model <- try(nls(coverage ~ L / (1 + exp(-k * (percentage - x0))),
                            data = agg_data,
                            start = list(L = 100, k = 0.1, x0 = 50),
                            control = nls.control(maxiter = 200)),
                        silent = TRUE)
  
  if (!inherits(logistic_model, "try-error")) {
    models$logistic <- logistic_model
  }
  
  # If all models failed, try a simple loess smoother
  if (length(models) == 0) {
    message("All nonlinear models failed to converge. Using loess smoothing instead.")
    loess_model <- loess(coverage ~ percentage, data = agg_data, span = 0.75)
    models$loess <- loess_model
  }
  
  # Create prediction grid for plotting
  pred_grid <- data.frame(percentage = seq(0, max(agg_data$percentage) * 1.5, length.out = 100))
  
  # Make predictions for each model
  for (model_name in names(models)) {
    if (model_name == "loess") {
      pred_grid[[model_name]] <- predict(models[[model_name]], newdata = pred_grid)
    } else {
      pred_grid[[model_name]] <- predict(models[[model_name]], newdata = pred_grid)
    }
  }
  
  # Create a plot with all models
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    library(reshape2)
    
    # Melt the prediction data for plotting
    pred_long <- melt(pred_grid, id.vars = "percentage", 
                      variable.name = "model", value.name = "predicted_coverage")
    
    # Get title based on coverage column used
    plot_title <- paste("Model Fitting for", gsub("_", " ", coverage_column))
    y_label <- paste("Percentage of", gsub("_", " ", coverage_column))
    
    # Create the plot
    p <- ggplot() +
      geom_point(data = agg_data, aes(x = percentage, y = coverage), alpha = 0.6, size = 3) +
      geom_line(data = pred_long, aes(x = percentage, y = predicted_coverage, color = model), size = 1) +
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
      scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 10))
    
    print(p)
  }
  
  # Estimate required percentage for target coverage
  result_table <- data.frame(
    model = character(),
    required_percentage = numeric(),
    max_achievable = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (model_name in names(models)) {
    # Predict coverage for a wide range of percentages
    extended_grid <- data.frame(percentage = seq(0, 1000, by = 0.1))
    
    if (model_name == "loess") {
      # For loess, we can only predict within the original data range
      extended_grid <- data.frame(percentage = seq(0, max(agg_data$percentage), by = 0.1))
      predicted <- predict(models[[model_name]], newdata = extended_grid)
    } else {
      predicted <- predict(models[[model_name]], newdata = extended_grid)
    }
    
    # Find max achievable coverage
    max_achievable <- max(predicted, na.rm = TRUE)
    
    # Calculate required percentage
    if (max_achievable >= target_coverage) {
      required_pct <- min(extended_grid$percentage[predicted >= target_coverage], na.rm = TRUE)
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
  
  # Print results
  print(result_table)
  
  return(list(
    models = models, 
    predictions = pred_grid, 
    required = result_table, 
    coverage_column = coverage_column
  ))
}
