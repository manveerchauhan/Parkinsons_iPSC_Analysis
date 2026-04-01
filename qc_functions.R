# QC Functions for Comprehensive Seurat Analysis
# Author: Manveer Chauhan
# These functions implement the full QC pipeline from FLAMES-seuratQC-script1.R

## Helper Function: Convert Cell Cycle Genes to Seurat Feature Format
convert_cc_genes_to_seurat_format <- function(gene_symbols, gene_dict_path) {

  # Load gene dictionary
  gene_dict <- read.csv(gene_dict_path, stringsAsFactors = FALSE)

  # Create ENSG-SYMBOL format from dictionary
  gene_dict <- gene_dict %>%
    mutate(feature_name = paste0(gene_id, "-", gene_symbol))

  # Match input gene symbols to dictionary
  matched_features <- gene_dict %>%
    filter(gene_symbol %in% gene_symbols) %>%
    pull(feature_name)

  # Report matching statistics
  n_input <- length(gene_symbols)
  n_matched <- length(matched_features)
  n_missing <- n_input - n_matched

  cat("  Gene matching: ", n_matched, "/", n_input, " genes matched (",
      round(100*n_matched/n_input, 1), "%)\n", sep="")

  if(n_missing > 0) {
    missing_genes <- setdiff(gene_symbols, gene_dict$gene_symbol)
    cat("  Missing genes:", paste(head(missing_genes, 10), collapse=", "))
    if(n_missing > 10) cat(" ... and", n_missing - 10, "more")
    cat("\n")
  }

  # Check if we have enough genes
  if(n_matched < n_input * 0.5) {
    warning("Less than 50% of cell cycle genes matched!")
  }

  return(matched_features)
}

## Function 2: Normalization and Variable Features
normalize_and_find_features <- function(seurat_obj, sample_id) {

  cat("\n=== STEP 2: Normalization and variable features for", sample_id, "===\n")

  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj)
  cat("  Normalization complete\n")

  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                     selection.method = 'vst',
                                     nfeatures = 2000)
  cat("  Found 2000 variable features\n")

  # Get top 10 variable features
  top10 <- head(VariableFeatures(seurat_obj), 10)

  # Create variable feature plot
  var_feature_plot <- VariableFeaturePlot(seurat_obj) +
    ggtitle(paste0(sample_id, ": Top 2000 Variable Genes")) +
    NoLegend()
  var_feature_plot <- LabelPoints(plot = var_feature_plot,
                                  points = top10,
                                  labels = top10,
                                  repel = TRUE)

  return(list(
    seurat_obj = seurat_obj,
    plots = list(var_feature_plot = var_feature_plot),
    top10 = top10
  ))
}

## Function 3: Cell Cycle Scoring and Scaling
cell_cycle_and_scale <- function(seurat_obj, sample_id, gene_dict_path = "./output_files/ref_files/isoform_gene_dict.csv") {

  cat("\n=== STEP 3: Cell cycle scoring and scaling for", sample_id, "===\n")

  # Convert cell cycle genes to ENSG-SYMBOL format
  cat("  Converting S-phase genes...\n")
  s_genes_converted <- convert_cc_genes_to_seurat_format(s_genes, gene_dict_path)

  cat("  Converting G2M-phase genes...\n")
  g2m_genes_converted <- convert_cc_genes_to_seurat_format(g2m_genes, gene_dict_path)

  # Cell cycle scoring with converted gene lists
  seurat_obj <- CellCycleScoring(seurat_obj,
                                 g2m.features = g2m_genes_converted,
                                 s.features = s_genes_converted)
  cat("  Cell cycle scoring complete\n")

  # Scale data (NO regression)
  all_features <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all_features)
  cat("  Data scaling complete\n")

  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  cat("  PCA complete\n")

  # Create cell cycle plots
  cc_pc1_pc2 <- DimPlot(seurat_obj, reduction = "pca",
                        label = FALSE, group.by = 'Phase') +
    ggtitle(paste0(sample_id, ": Cell Cycle on PC1 & PC2")) +
    labs(color = "Cell Phase")

  cc_pc2_pc3 <- DimPlot(seurat_obj, reduction = "pca",
                        dims = c(2,3),
                        label = FALSE, group.by = 'Phase') +
    ggtitle(paste0(sample_id, ": Cell Cycle on PC2 & PC3")) +
    labs(color = "Cell Phase")

  return(list(
    seurat_obj = seurat_obj,
    plots = list(cc_pc1_pc2 = cc_pc1_pc2, cc_pc2_pc3 = cc_pc2_pc3)
  ))
}

## Function 4: Quantitative Elbow Analysis
quantitative_elbow <- function(seurat_obj, sample_id, max_pcs = 50, pc_adjustment = 5) {

  cat("\n=== STEP 4: Quantitative elbow analysis for", sample_id, "===\n")

  # Calculate percent variance explained by each PC
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)

  # Metric 1: PC where cumulative >90% AND individual <5%
  co1 <- which(cumu > 90 & pct < 5)[1]

  # Metric 2: PC where consecutive change <0.1%
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),
              decreasing = TRUE)[1] + 1

  # Optimal PCs (calculated)
  optimal_pcs_calculated <- min(co1, co2, na.rm = TRUE)

  # Adjusted PCs (calculated + adjustment)
  optimal_pcs_used <- optimal_pcs_calculated + pc_adjustment

  cat("  Metric 1 (cumu >90% & pct <5%):", co1, "PCs\n")
  cat("  Metric 2 (consecutive change <0.1%):", co2, "PCs\n")
  cat("  Optimal PCs calculated:", optimal_pcs_calculated, "\n")
  cat("  Optimal PCs used (calculated +", pc_adjustment, "):", optimal_pcs_used, "\n")

  # Create elbow plot with subtitle showing both values
  elbow_plot <- ElbowPlot(seurat_obj, ndims = max_pcs) +
    ggtitle(paste0(sample_id, ": Elbow Plot")) +
    labs(subtitle = paste0("Calculated: ", optimal_pcs_calculated, " PCs | ",
                          "Used: ", optimal_pcs_used, " PCs (+", pc_adjustment, ")")) +
    geom_vline(xintercept = optimal_pcs_calculated, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_vline(xintercept = optimal_pcs_used, linetype = "solid", color = "darkgreen", linewidth = 1) +
    annotate("text", x = optimal_pcs_calculated, y = Inf,
             label = paste0("Calculated: ", optimal_pcs_calculated),
             vjust = 1.5, hjust = -0.1, size = 3, color = "red") +
    annotate("text", x = optimal_pcs_used, y = Inf,
             label = paste0("Used: ", optimal_pcs_used),
             vjust = 3, hjust = -0.1, size = 3, color = "darkgreen")

  # Create quantitative plot with subtitle
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  quant_plot <- ggplot(plot_df, aes(cumu, pct, label = rank,
                                    color = rank > optimal_pcs_calculated)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    ggtitle(paste0(sample_id, ": Quantitative PC Selection")) +
    labs(subtitle = paste0("Calculated: ", optimal_pcs_calculated, " PCs | ",
                          "Used: ", optimal_pcs_used, " PCs (+", pc_adjustment, ")")) +
    xlab("Cumulative Variance (%)") +
    ylab("Variance (%)") +
    theme_bw()

  return(list(
    optimal_pcs = optimal_pcs_calculated,  # Return calculated value for backward compatibility
    optimal_pcs_used = optimal_pcs_used,   # Return adjusted value
    plots = list(elbow_plot = elbow_plot, quant_plot = quant_plot),
    metrics = list(co1 = co1, co2 = co2, pct = pct, cumu = cumu)
  ))
}

## Function 5: Clustering Optimization with Silhouette Analysis
optimize_clustering <- function(seurat_obj, sample_id, optimal_pcs) {

  cat("\n=== STEP 5: Clustering optimization for", sample_id, "===\n")

  # Find neighbors
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:optimal_pcs)
  cat("  FindNeighbors complete (dims 1:", optimal_pcs, ")\n")

  # Silhouette analysis
  cat("  Running silhouette analysis...\n")
  sil_results <- optimize_silhouette(sobject = seurat_obj,
                                     test_res = seq(0.1, 1.2, by = 0.1),
                                     summary_plot = TRUE,
                                     reduction = "pca")

  sil_df <- as.data.frame(sil_results) %>%
    dplyr::rename(avg_sil_vals = sil_vals) %>%
    group_by(num_clusters) %>%
    dplyr::slice(which.max(avg_sil_vals)) %>%
    ungroup() %>%
    arrange(desc(avg_sil_vals))

  # Get optimal resolution
  optimal_res <- sil_df$res_vals[1]
  cat("  Optimal resolution:", optimal_res, "(silhouette:",
      round(sil_df$avg_sil_vals[1], 4), ")\n")

  # Cluster with all test resolutions for clustree
  seurat_obj <- FindClusters(seurat_obj, resolution = seq(0.1, 1.2, by = 0.1))

  # Create clustree plot
  clustree_plot <- clustree(seurat_obj) +
    ggtitle(paste0(sample_id, ": Clustree"))

  # Re-cluster with optimal resolution
  seurat_obj <- FindClusters(seurat_obj, resolution = optimal_res)

  return(list(
    seurat_obj = seurat_obj,
    optimal_res = optimal_res,
    sil_results = sil_df,
    plots = list(clustree = clustree_plot)
  ))
}

## Function 6: UMAP and Doublet Detection
umap_and_doublets <- function(seurat_obj, sample_id, optimal_pcs, doublet_rate) {

  cat("\n=== STEP 6: UMAP and doublet detection for", sample_id, "===\n")

  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:optimal_pcs)
  cat("  UMAP complete\n")

  # DoubletFinder parameter sweep
  cat("  Running DoubletFinder parameter sweep...\n")
  sweep_res_list <- paramSweep(seurat_obj, PCs = 1:optimal_pcs, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  BCmvn <- find.pK(sweep_stats)

  # Select optimal pK
  pK <- BCmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  cat("  Optimal pK:", pK, "\n")

  # Homotypic doublet proportion
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)

  # Calculate expected doublets
  nExp_poi <- round(doublet_rate * nrow(seurat_obj@meta.data))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  cat("  Expected doublets:", nExp_poi_adj, "\n")

  # Run DoubletFinder
  seurat_obj <- doubletFinder(seurat_obj,
                              PCs = 1:optimal_pcs,
                              pN = 0.25,
                              pK = pK,
                              nExp = nExp_poi_adj,
                              reuse.pANN = FALSE,
                              sct = FALSE)

  # Clean up column names
  colnames(seurat_obj@meta.data) <- sub("DF.classifications_.*$",
                                         "DF.classifications",
                                         colnames(seurat_obj@meta.data))

  # Doublet statistics
  doublet_stats <- seurat_obj@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(median_nCount = median(nCount_RNA),
              median_nFeature = median(nFeature_RNA),
              count = n())

  cat("  Doublets detected:\n")
  print(doublet_stats)

  # Plot doublets
  doublet_umap <- DimPlot(seurat_obj, reduction = 'umap', group.by = "DF.classifications") +
    ggtitle(paste0(sample_id, ": Doublets Detected"))

  # Remove doublets
  cells_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = DF.classifications == 'Singlet')
  cells_after <- ncol(seurat_obj)
  cat("  Cells before doublet removal:", cells_before, "\n")
  cat("  Cells after doublet removal:", cells_after, "\n")

  # Final UMAP
  final_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle(paste0(sample_id, ": Final UMAP")) +
    labs(color = "Cluster")

  # QC metrics UMAP plots
  nfeature_umap <- FeaturePlot(seurat_obj, reduction = "umap", features = "nFeature_RNA") +
    ggtitle(paste0(sample_id, ": nFeature_RNA"))

  ncount_umap <- FeaturePlot(seurat_obj, reduction = "umap", features = "nCount_RNA") +
    ggtitle(paste0(sample_id, ": nCount_RNA"))

  # Cell cycle UMAP plot
  cellcycle_umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "Phase") +
    ggtitle(paste0(sample_id, ": Cell Cycle Phase")) +
    labs(color = "Phase")

  return(list(
    seurat_obj = seurat_obj,
    doublet_stats = doublet_stats,
    plots = list(
      doublet_umap = doublet_umap,
      final_umap = final_umap,
      nfeature_umap = nfeature_umap,
      ncount_umap = ncount_umap,
      cellcycle_umap = cellcycle_umap
    )
  ))
}
