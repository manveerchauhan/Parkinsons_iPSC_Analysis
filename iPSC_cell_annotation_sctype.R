## Automated Cell Type Annotation using scType------
# Author: Manveer Chauhan
# This script performs automated cell type annotation on the Harmony-integrated
# iPSC Parkinson's object using scType with ScTypeDB_Default and PanglaoDB references.
#
# Adapted from: Neurodevelopmental_Models_analysis/draft_space/ch3_cell_annotation_sctype.R
# Input  : output_files/integrated_objects/iPSC_PD_integrated_harmony.rds
# Output : output_files/cell_annotation/

# Load essential packages
library(Seurat)
library(tidyverse)
library(grid)
library(gridExtra)
library(openxlsx)
library(HGNChelper)
library(patchwork)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Parkinsons_iPSC_Analysis")

# Create output directory
dir.create("./output_files/cell_annotation", recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(5728)

cat("\n========================================\n")
cat("SCTYPE CELL ANNOTATION PIPELINE — iPSC PD\n")
cat("========================================\n\n")

# Load scType's gene set preparation function
cat("Loading scType functions from GitHub...\n")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# Load scType's cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
cat("scType functions loaded successfully\n\n")

## DEFINE REFERENCE DATABASES------
cat("=== Configuring reference databases ===\n")

reference_databases <- list(
  "ScTypeDB_Default" = list(
    use_alt = FALSE,
    filepath = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
    description = "Default ScTypeDB Brain markers"
  ),
  "Panglao" = list(
    use_alt = TRUE,
    filepath = "/data/gpfs/projects/punim0646/manveer/panglao-sctypeFormat.xlsx",
    description = "PanglaoDB markers (scType format)"
  )
)

cat("Configured", length(reference_databases), "reference databases:\n")
for(ref_name in names(reference_databases)) {
  cat("  -", ref_name, ":", reference_databases[[ref_name]]$description, "\n")
}
cat("\n")

## FUNCTION 1: Gene Name Resolution (Reverse Mapping)------
resolve_sctype_genes <- function(sctype_gene_list, seurat_obj) {
  # Maps scType gene symbols to ENSG-SYMBOL format features in Seurat object
  # This is the reverse mapping approach - we map scType genes to our format

  cat("    Resolving scType gene symbols to ENSG-SYMBOL format...\n")

  # Get all features from Seurat object
  all_features <- rownames(seurat_obj)

  # Extract symbols: ENSG00000123456.10_CD99 -> CD99
  # Pattern matches everything after the last underscore or dash
  feature_symbols <- sub(".*[-_]", "", all_features)

  # Create mapping: gene_symbol -> full_feature_name
  symbol_to_feature <- setNames(all_features, feature_symbols)

  # Map scType's positive gene lists to full feature names
  resolved_positive <- list()
  total_positive_genes <- 0
  total_positive_resolved <- 0

  for(cell_type in names(sctype_gene_list$gs_positive)) {
    original_genes <- sctype_gene_list$gs_positive[[cell_type]]
    resolved <- symbol_to_feature[original_genes]
    resolved <- resolved[!is.na(resolved)]  # Remove unmapped genes
    resolved_positive[[cell_type]] <- unname(resolved)

    total_positive_genes <- total_positive_genes + length(original_genes)
    total_positive_resolved <- total_positive_resolved + length(resolved)
  }

  # Map scType's negative gene lists to full feature names
  resolved_negative <- list()
  total_negative_genes <- 0
  total_negative_resolved <- 0

  for(cell_type in names(sctype_gene_list$gs_negative)) {
    original_genes <- sctype_gene_list$gs_negative[[cell_type]]
    resolved <- symbol_to_feature[original_genes]
    resolved <- resolved[!is.na(resolved)]
    resolved_negative[[cell_type]] <- unname(resolved)

    total_negative_genes <- total_negative_genes + length(original_genes)
    total_negative_resolved <- total_negative_resolved + length(resolved)
  }

  # Report statistics
  pct_positive <- round(100 * total_positive_resolved / total_positive_genes, 1)
  pct_negative <- round(100 * total_negative_resolved / total_negative_genes, 1)

  cat("    Positive markers:", total_positive_resolved, "/", total_positive_genes,
      "(", pct_positive, "%)\n")
  cat("    Negative markers:", total_negative_resolved, "/", total_negative_genes,
      "(", pct_negative, "%)\n")

  return(list(gs_positive = resolved_positive, gs_negative = resolved_negative))
}

## FUNCTION 2: Main Annotation Function------
annotate_cellTypes_sctype <- function(seurat_obj, timepoint_name, reference_name,
                                      reference_info) {

  cat("\n  === Annotating with", reference_name, "===\n")
  cat("  Description:", reference_info$description, "\n")

  # Load the marker database
  db_ <- reference_info$filepath
  tissue <- "Brain"

  cat("  Loading marker database...\n")

  # Prepare gene sets
  gs_list <- gene_sets_prepare(db_, tissue)
  cat("  Loaded", length(gs_list$gs_positive), "cell type markers\n")

  # Resolve gene symbols to ENSG-SYMBOL format
  gs_list_resolved <- resolve_sctype_genes(gs_list, seurat_obj)

  # Get cell-type by cell matrix using resolved gene lists
  cat("  Calculating scType scores...\n")
  es.max <- sctype_score(scRNAseqData = seurat_obj[["RNA"]]$scale.data,
                         scaled = TRUE,
                         gs = gs_list_resolved$gs_positive,
                         gs2 = gs_list_resolved$gs_negative)

  # Merge by cluster
  cat("  Assigning cell types to clusters...\n")
  cL_results <- do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
    es.max.cl <- sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]),
                      decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl,
                    ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
  }))

  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)

  # Set low-confidence clusters to "Unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

  cat("  Cluster annotations:\n")
  print(sctype_scores[,1:3])

  # Add annotations to metadata
  annotation_col <- paste0("sctype_", reference_name)
  seurat_obj@meta.data[[annotation_col]] <- ""

  for(j in unique(sctype_scores$cluster)){
    cl_type <- sctype_scores[sctype_scores$cluster==j,]
    seurat_obj@meta.data[[annotation_col]][seurat_obj@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
  }

  # Create UMAPs
  cat("  Creating UMAP plots...\n")

  labeledUMAP <- DimPlot(seurat_obj,
                         reduction = "umap.harmony",
                         label = TRUE,
                         label.size = 3,
                         repel = TRUE,
                         group.by = annotation_col) +
    labs(title = paste0(timepoint_name, ": Cell-type Annotations"),
         subtitle = paste0("Reference: ", reference_name),
         color = "Cell Type") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10))

  baseUMAP <- DimPlot(seurat_obj,
                      reduction = "umap.harmony",
                      label = TRUE,
                      repel = TRUE) +
    labs(title = paste0(timepoint_name, ": Unsupervised Clustering"),
         subtitle = "From Harmony integration",
         color = "Cluster") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10))

  # Create pie charts for confidence scores
  cat("  Creating confidence pie charts...\n")

  # Prepare the edges
  cL_results <- cL_results[order(cL_results$cluster),]
  edges <- cL_results
  edges$type <- paste0(edges$type, "_", edges$cluster)
  edges$cluster <- paste0("cluster ", edges$cluster)
  edges <- edges[, c("cluster", "type")]
  colnames(edges) <- c("from", "to")
  rownames(edges) <- NULL

  # Prepare nodes
  nodes_lvl1 <- sctype_scores[, c("cluster", "ncells")]
  nodes_lvl1$cluster <- paste0("cluster ", nodes_lvl1$cluster)
  nodes_lvl1$Colour <- "#f1f1ef"
  nodes_lvl1$ord <- 1
  nodes_lvl1$realname <- nodes_lvl1$cluster
  nodes_lvl1 <- as.data.frame(nodes_lvl1)

  nodes_lvl2 <- data.frame(cluster = character(),
                           ncells = numeric(),
                           Colour = character(),
                           ord = numeric(),
                           realname = character(),
                           stringsAsFactors = FALSE)

  ccolss <- c("#5f75ae", "#92bbb8", "#64a841", "#e5486e", "#de8e06", "#eccf5a",
              "#b5aa0f", "#e4b680", "#7ba39d", "#b15928", "#ffff99", "#6a3d9a",
              "#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99", "#33a02c",
              "#b2df8a", "#1f78b4", "#a6cee3")

  for (i in 1:length(unique(cL_results$cluster))) {
    dt_tmp <- cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]
    nodes_tmp <- data.frame(cluster = paste0(dt_tmp$type, "_", dt_tmp$cluster),
                            ncells = dt_tmp$scores,
                            Colour = ccolss[i],
                            ord = 2,
                            realname = dt_tmp$type,
                            stringsAsFactors = FALSE)
    nodes_lvl2 <- rbind(nodes_lvl2, nodes_tmp)
  }

  nodes <- rbind(nodes_lvl1, nodes_lvl2)
  nodes$ncells[nodes$ncells < 1] <- 1

  files_db <- openxlsx::read.xlsx(db_)[, c("cellName", "shortName")]
  files_db <- unique(files_db)
  nodes <- merge(nodes, files_db, all.x = TRUE, all.y = FALSE,
                 by.x = "realname", by.y = "cellName", sort = FALSE)
  nodes$shortName[is.na(nodes$shortName)] <- nodes$realname[is.na(nodes$shortName)]
  nodes <- nodes[, c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

  # Remove duplicates
  duplicated_rows <- duplicated(nodes$cluster)
  nodes <- nodes[!duplicated_rows, ]

  # Filter for pie charts
  filtered_nodes <- subset(nodes, ord != 1)
  filtered_nodes <- subset(filtered_nodes, select = -c(4, 5))

  # Extract cluster numbers
  filtered_nodes$cluster <- str_replace(filtered_nodes$cluster, ".*_([0-9]+)", "Cluster \\1")

  # Calculate proportions
  grouped_df <- filtered_nodes %>%
    group_by(cluster) %>%
    mutate(prop = ncells / sum(ncells)) %>%
    ungroup() %>%
    arrange(cluster)

  # Create pie charts
  pieChartList <- list()
  for (cluster in unique(grouped_df$cluster)) {
    cluster_data <- grouped_df[grouped_df$cluster == cluster, ]

    pie_chart <- ggplot(cluster_data, aes(x = "", y = prop, fill = realname)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      labs(title = paste(cluster),
           subtitle = paste0("Reference: ", reference_name)) +
      theme_void() +
      theme(legend.position = "right",
            plot.subtitle = element_text(size = 8)) +
      labs(fill = "Cell Types")

    pie_chart <- pie_chart +
      geom_text(aes(label = paste0(round(prop * 100, 1), "%"), x = 1.5),
                position = position_stack(vjust = 0.5), size = 3)

    pieChartList <- c(pieChartList, list(pie_chart))
  }

  cat("  Created", length(pieChartList), "pie charts\n")

  # Return results
  return(list(
    seurat_obj = seurat_obj,
    baseUMAP = baseUMAP,
    labeledUMAP = labeledUMAP,
    pieChartList = pieChartList,
    sctype_scores = sctype_scores,
    annotation_col = annotation_col
  ))
}

## FUNCTION 3: Helper Function for Pie Chart Layout------
closest_factor_pair <- function(n) {
  # Find the closest factor pair for grid layout
  sqrt_floor <- floor(sqrt(n))

  closest_factor1 <- 1
  closest_factor2 <- n

  for (i in sqrt_floor:1) {
    if (n %% i == 0) {
      closest_factor1 <- i
      closest_factor2 <- n / i
      break
    }
  }

  factor.pair <- c(closest_factor1, closest_factor2)
  return(factor.pair)
}

## FUNCTION 4: Generate Pie Chart PDF------
generate_pieChartPDFs <- function(annotation_result, timepoint_name, reference_name) {

  piechartlist <- annotation_result$pieChartList

  cluster.num <- length(piechartlist)
  row.num <- closest_factor_pair(cluster.num)[2]
  col.num <- closest_factor_pair(cluster.num)[1]
  pdf.width <- 7*col.num
  pdf.height <- 5*row.num

  # Create filename
  fileName <- paste0("./output_files/cell_annotation/",
                     timepoint_name, "_", reference_name, "_confidence_piecharts.pdf")

  cat("  Saving pie charts:", fileName, "\n")

  pdf(fileName, width = pdf.width, height = pdf.height)

  # Add title
  title_text <- paste0(timepoint_name, " - scType Confidence Scores\n",
                       "Reference: ", reference_name)

  piechart.layout <- grid.arrange(grobs = piechartlist,
                                   nrow = row.num,
                                   ncol = col.num,
                                   top = textGrob(title_text,
                                                  gp = gpar(fontsize = 16, fontface = "bold")))

  dev.off()

  return(fileName)
}

## FUNCTION 5: Generate UMAP PDF------
generate_UMAP_PDFs <- function(annotation_result, timepoint_name, reference_name) {

  fileName <- paste0("./output_files/cell_annotation/",
                     timepoint_name, "_", reference_name, "_UMAPs.pdf")

  cat("  Saving UMAPs:", fileName, "\n")

  pdf(fileName, width = 14, height = 6)

  # Create side-by-side layout
  combined_plot <- annotation_result$labeledUMAP | annotation_result$baseUMAP
  print(combined_plot)

  dev.off()

  return(fileName)
}

## FUNCTION 6: Create Annotation Summary Table------
create_annotation_summary_table <- function(all_results, seurat_obj, timepoint_name) {

  cat("\n  === Creating annotation summary table ===\n")

  # Get unique clusters
  clusters <- sort(unique(seurat_obj$seurat_clusters))

  # Initialize summary dataframe
  summary_df <- data.frame(
    Timepoint = timepoint_name,
    Cluster = clusters,
    stringsAsFactors = FALSE
  )

  # Add cell counts
  cluster_counts <- table(seurat_obj$seurat_clusters)
  summary_df$nCells <- as.numeric(cluster_counts[as.character(summary_df$Cluster)])

  # Add annotation from each reference
  for(ref_name in names(all_results)) {
    annotations <- all_results[[ref_name]]$sctype_scores

    # Match cluster annotations
    for(i in 1:nrow(summary_df)) {
      cluster_id <- summary_df$Cluster[i]
      annotation_row <- annotations[annotations$cluster == cluster_id, ]

      if(nrow(annotation_row) > 0) {
        summary_df[[ref_name]][i] <- as.character(annotation_row$type[1])
      } else {
        summary_df[[ref_name]][i] <- "Unknown"
      }
    }
  }

  cat("  Summary table created with", nrow(summary_df), "clusters\n")

  return(summary_df)
}

## FUNCTION 7: Generate Summary Report PDF------
generate_summary_report_pdf <- function(summary_df, all_results, timepoint_name) {

  cat("\n  === Generating summary report PDF ===\n")

  fileName <- paste0("./output_files/cell_annotation/",
                     timepoint_name, "_annotation_summary.pdf")

  pdf(fileName, width = 16, height = 10)

  # Page 1: Summary Table
  title_text <- paste0(timepoint_name, " - Annotation Summary Across All References")

  # Create table grob
  table_grob <- tableGrob(summary_df, rows = NULL,
                          theme = ttheme_default(base_size = 9))

  grid.arrange(table_grob,
               top = textGrob(title_text,
                             gp = gpar(fontsize = 14, fontface = "bold")))

  # Page 2: Side-by-side UMAP comparison (2x2 grid)
  umap_list <- list()
  for(ref_name in names(all_results)) {
    umap_list[[ref_name]] <- all_results[[ref_name]]$labeledUMAP +
      theme(legend.position = "none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
  }

  # Create grid of reference UMAPs
  comparison_title <- paste0(timepoint_name, " - UMAP Comparison Across References")

  combined_umaps <- wrap_plots(umap_list, ncol = 2) +
    plot_annotation(title = comparison_title,
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  print(combined_umaps)

  dev.off()

  cat("  Summary report saved:", fileName, "\n")

  return(fileName)
}

## FUNCTION 8: Process One Timepoint with All References------
process_timepoint_annotation <- function(timepoint_name, seurat_obj, reference_list) {

  cat("\n\n####################################################\n")
  cat("PROCESSING:", timepoint_name, "\n")
  cat("####################################################\n\n")

  cat("Loaded object:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "features\n")
  cat("Number of clusters:", length(unique(seurat_obj$seurat_clusters)), "\n\n")

  # Join layers for Seurat v5 compatibility
  cat("Joining layers for Seurat v5 compatibility...\n")
  seurat_obj <- JoinLayers(seurat_obj)
  cat("Layers joined successfully\n")

  # Process with each reference database
  all_results <- list()

  for(ref_name in names(reference_list)) {

    ref_info <- reference_list[[ref_name]]

    # Run annotation
    annotation_result <- annotate_cellTypes_sctype(
      seurat_obj = seurat_obj,
      timepoint_name = timepoint_name,
      reference_name = ref_name,
      reference_info = ref_info
    )

    # Update seurat object with new annotations
    seurat_obj <- annotation_result$seurat_obj

    # Generate PDFs
    generate_UMAP_PDFs(annotation_result, timepoint_name, ref_name)
    generate_pieChartPDFs(annotation_result, timepoint_name, ref_name)

    # Store results
    all_results[[ref_name]] <- annotation_result

    cat("\n")
  }

  # Create summary table
  summary_df <- create_annotation_summary_table(all_results, seurat_obj, timepoint_name)

  # Save summary CSV
  csv_file <- paste0("./output_files/cell_annotation/",
                     timepoint_name, "_annotation_summary.csv")
  write.csv(summary_df, csv_file, row.names = FALSE)
  cat("\n  Summary CSV saved:", csv_file, "\n")

  # Generate summary report PDF
  generate_summary_report_pdf(summary_df, all_results, timepoint_name)

  cat("\n====================================\n")
  cat("ANNOTATION COMPLETE:", timepoint_name, "\n")
  cat("====================================\n")

  return(list(
    seurat_obj = seurat_obj,
    all_results = all_results,
    summary_df = summary_df
  ))
}

## MAIN EXECUTION------
cat("\n\n####################################################\n")
cat("STARTING CELL ANNOTATION — iPSC PD\n")
cat("####################################################\n\n")

INPUT_RDS  <- "output_files/integrated_objects/iPSC_PD_integrated_harmony.rds"
OUTPUT_RDS <- "output_files/integrated_objects/iPSC_PD_integrated_harmony_annotated.rds"

# Load integrated object
cat("Loading integrated object:", INPUT_RDS, "\n")
if (!file.exists(INPUT_RDS)) {
  stop("Could not find: ", INPUT_RDS, "\n  Run iPSC_harmony_integration.R first.")
}
seurat_obj <- readRDS(INPUT_RDS)

# Run annotation with all references
results <- process_timepoint_annotation(
  timepoint_name = "iPSC_PD",
  seurat_obj = seurat_obj,
  reference_list = reference_databases
)

# Save annotated object
saveRDS(results$seurat_obj, file = OUTPUT_RDS)
cat("\nAnnotated object saved:", OUTPUT_RDS, "\n")

# Force garbage collection
gc()

cat("\n\n####################################################\n")
cat("ANNOTATION COMPLETE — iPSC PD\n")
cat("####################################################\n")

cat("\n=== ANNOTATION SUMMARY ===\n")
cat("  Total cells     :", ncol(results$seurat_obj), "\n")
cat("  Number of clusters:", length(unique(results$seurat_obj$seurat_clusters)), "\n")
cat("  References tested:", length(reference_databases), "\n")
cat("  Metadata columns added:", length(reference_databases), "\n")

cat("\n=== ALL ANNOTATION COMPLETE ===\n")
cat("Annotated object saved in: output_files/integrated_objects/\n")
cat("PDF reports and CSVs saved in: output_files/cell_annotation/\n")
