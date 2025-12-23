library(DESeq2)
library(ggplot2)
library(glue)
library(EnhancedVolcano)
library(org.Ce.eg.db)
library(ComplexHeatmap)
library(GeneOverlap)
library(ggforce)
library(factoextra)
library(extrafont)
library(ggpubr)
library(tibble)
library(ReportingTools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dbscan)
library(waRRior)
library(clusterProfiler)
library(dplyr)
#font_import()
#loadfonts(device = "win")
data(GeneOverlap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~constants~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
readsLim <- 5
lockBinding("readsLim", globalenv())

alphaLim <- 0.05
lockBinding("alphaLim", globalenv())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prep_counts <- function(file) {
  counts <- read.delim(file, header=TRUE, row.names = 1, sep = ",")
  counts <- counts[which(rowSums(counts) >= readsLim), ]
  return(counts)
}

prep_coldata <- function(file) {
  # Set stringsAsFactors to TRUE for automatic factor conversion
  design_file <- read.csv(file, sep=',', row.names = 1, stringsAsFactors = TRUE)
  return(design_file)
}

order_data <- function(coldata, df) {
  # Define the desired order of conditions
  desired_order <- c("start", "intermediate", "end")
  # Extract row names based on conditions
  start_rows <- rownames(coldata)[coldata$condition == "start"]
  intermediate_rows <- rownames(coldata)[coldata$condition == "intermediate"]
  end_rows <- rownames(coldata)[coldata$condition == "end"]
  # Order the data frame columns based on the desired order
  desired_cols <- c(start_rows, intermediate_rows, end_rows)
  ordered_df <- df[, desired_cols]
  return(ordered_df)
}

prepare_data <- function(counts, refer, coldata) {
  ordered <- order_data(coldata, counts)
  ordered_coldata <- coldata[colnames(ordered), ]
  
  # Create the DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = ordered, colData = ordered_coldata, design = ~condition + worm)
  dds$condition <- relevel(dds$condition, ref = refer)
  
  # Run DESeq analysis with local fit to avoid warnings about parametric fit
  dds <- DESeq(dds)
  
  return(dds)
}

# Perform PCA (Principal Component Analysis) on rlog-transformed data. 
# A modified versions of the opensource code of built-in deseq2 functions for better visualization
do_pca <- function(dds, pc1, pc2, intgroup="condition", p_colors = c("blue","black", "seagreen"), ntop=500, returnData=FALSE) {
  transformed <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5))
  
  rv <- rowVars(assay(transformed))  # Compute row variances
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]  # Select top ntop variable genes
  
  pca <- prcomp(t(assay(transformed)[select,]))  # Perform PCA on selected genes
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)  # Compute variance explained by each principal component
  
  # Ensure that the columns specified by intgroup exist in the colData
  if (!all(intgroup %in% names(colData(transformed)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(transformed)[, intgroup, drop=FALSE])  # Extract colData for grouping
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse=":"))  # Concatenate multiple groups into one factor
  } else {
    colData(transformed)[[intgroup]]  # Use single grouping variable
  }
  
  # Create a data frame for plotting PCA results
  d <- data.frame(PC1=pca$x[,pc1], PC2=pca$x[,pc2], group=group, intgroup.df, name=colnames(transformed))
  
  if (returnData) {  # If returnData is TRUE, return the data frame instead of plotting
    attr(d, "percentVar") <- percentVar[c(pc1, pc2)]
    return(d)
  }
  
  # Create PCA plot using ggplot2
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
    geom_point(size=3) +  # Plot points for each sample
    geom_text(aes(label=name), size=3, vjust=-0.5) +  # Add sample names as labels directly to each point
    xlab(paste0("PC",pc1," (",round(percentVar[pc1] * 100),"%)")) +
    ylab(paste0("PC",pc2," (",round(percentVar[pc2] * 100),"%)")) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.title = element_blank(), text=element_text(size=12,  family="sans")) + 
    coord_fixed() +  # Removed xlim and ylim to allow automatic scaling
    scale_color_manual(values = p_colors) +  # Customize colors
    scale_fill_manual(values = p_colors)
}

#look at the results of dds
results_data <- function(dds, toCompare, baseLevel) {
  res <- results(dds, alpha = alphaLim, contrast = c("condition", toCompare, baseLevel)) 
  head(results(dds, tidy=TRUE))
  res <- res[order(res$padj),]
  head(res,10)
  summary(res)
  # write.csv(as.data.frame(res), file=glue('{toCompare}_vs_{baseLevel}_DGE.csv'))
  return(res)
}

results_stable_genes <- function(dds, toCompare, baseLevel, lfcT) {
  res <- results(dds, alpha = alphaLim, lfcThreshold = lfcT, altHypothesis = "lessAbs", contrast = c("condition", toCompare, baseLevel)) 
  res <- res[order(res$padj),]
  res <- na.omit(res)
  head(res,10)
  summary(res)
  # write.csv(as.data.frame(res), file=glue('{toCompare}_vs_{baseLevel}_stable_genes.csv'))
  return(res)
}

do_heatmap_means <- function(res, dds, coldata, compared_condition, col = "slice") {
  # Filter significant genes
  sigs <- na.omit(res) 
  sigs <- sigs[sigs$padj < alphaLim, ]
  sigs.df <- as.data.frame(sigs) 
  sigs.df <- sigs.df[(abs(sigs.df$log2FoldChange) > 0) ,] 
  sigs.df <- sigs.df[order(sigs.df$log2FoldChange, decreasing = TRUE),]
  
  # Perform variance stabilizing transformation
  normalized_out <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5))
  
  # Filter for the relevant conditions
  selected_rows <- coldata[coldata$condition %in% compared_condition, ]
  
  # Extract the relevant genes and samples from the normalized matrix
  mat <- assay(normalized_out)[rownames(sigs.df), , drop = FALSE]
  
  # Initialize a list to store row means for each unique value in the specified column
  value_means <- list()
  
  # Loop through each unique value in the specified column
  for (val in unique(selected_rows[[col]])) {
    # Extract indices for the current value across all worms
    value_indices <- which(selected_rows[[col]] == val)
    
    # Combine worm, slice, and slice_index information to form the full column names
    full_column_names <- paste0(selected_rows$worm[value_indices], "_", selected_rows$slice[value_indices], "_", selected_rows$slice_index[value_indices])
    
    # Extract data for the current value
    value_data <- mat[, full_column_names, drop = FALSE]
    
    # Calculate row means for the current value
    value_mean <- rowMeans(value_data)
    
    # Assign a unique name to the column based on the value
    col_name <- paste0("Value_", val)
    
    # Store row means for the current value in the list
    value_means[[col_name]] <- value_mean
  }
  
  # Combine row means for all values into a single matrix
  mat_means <- do.call(cbind, value_means)
  
  # Center and scale each column (Z-score) then transpose
  mat.scaled <- t(apply(mat_means, 1, scale)) 
  
  mat.scaled <- pmax(pmin(mat.scaled, 2.5), -2.5)
  
  # Set column names for the scaled matrix
  colnames(mat.scaled) <- names(value_means)
  
  # Create the heatmap
  h1 <- Heatmap(mat.scaled, km=3, cluster_rows = F, cluster_columns = F, 
                column_labels = colnames(mat.scaled), use_raster = FALSE, name = "z_score", show_row_names = FALSE,
                heatmap_legend_param = list(at = seq(-2.5, 2.5, by = 1)))
  
  # Print the heatmap
  print(h1)
}

#creating heat map
# choosing significantly differential expressed genes, ordering them by expression,
# and than calculating the z score from the counts
do_heatmap <- function(res, dds, coldata, compared_condition) {
  
  sigs <- na.omit(res) # omitting NA from res
  sigs <- sigs[sigs$padj < alphaLim, ] # taking only significantly expressed genes
  sigs.df <- as.data.frame(sigs) 
  sigs.df <- sigs.df[(abs(sigs.df$log2FoldChange) > 0),] # arbitrary thresholds, can add & (sigs.df$baseMean > num)
  sigs.df <- sigs.df[order(sigs.df$log2FoldChange, decreasing = TRUE),] # order by log2foldchange
  
  normalized_out <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5)) # get normalized count data from dds object
  selected_rows <- coldata[coldata$condition %in% compared_condition, ] # choosing only the compared samples
  mat<-assay(normalized_out)[rownames(sigs.df), rownames(selected_rows)] # sig genes x samples 
  colnames(mat) <- rownames(selected_rows)
  mat.scaled <- t(apply(mat, 1, scale)) # center and scale each column (Z-score) then transpose
  colnames(mat.scaled)<-colnames(mat)
  
  
  #colHeat <- colorRamp2(c(-1, 0, 1), c("blue4","white", "red3")) and add col = colHeat
  
  h1 <- Heatmap(mat.scaled,km=3, use_raster = F, cluster_rows = F, cluster_columns = F, 
                column_labels = colnames(mat.scaled), name = "z_score", show_row_names = FALSE)
  
  #tiff("Heatmap_wt_FBF2_germline_cluster_by_samples_and genes.tiff", width = 1000, height = 1200, res = 150)
  print(h1)
  #dev.off()
}

#count up
count_upregulated <- function(res) {
  res <- na.omit(res)
  up <- res[(res$log2FoldChange>0)&(res$padj<0.05),]
  return(nrow(up))
}

#count down
count_downregulated <- function(res) {
  res <- na.omit(res)
  down <- res[(res$log2FoldChange<0)&(res$padj<0.05),]
  return(nrow(down))
}

#extracting upregulated genes
extract_upregulated_genes <- function(res, refer, compareTo) {
  res <- na.omit(res)
  up <- res[(res$log2FoldChange>0)&(res$padj<0.05),]
  #write.csv(as.data.frame(up), file=glue('Upregulated_{refer}_vs_{compareTo}.csv'))
  return(up)
}

#extracting downregulated genes
extract_downregulated_genes <- function(res, refer, compareTo) {
  res <- na.omit(res)
  down <- res[(res$log2FoldChange<0)&(res$padj<0.05),]
  #write.csv(as.data.frame(down), file=glue('Downregulated_{refer}_vs_{compareTo}.csv'))
  return(down)
}

myggmaplot <- function (data, subsetfile, subsetcolor1, subsetcolor2, fdr = 0.05, fc = 1.5, genenames = NULL,
                        detection_call = NULL, size = NULL, alpha = 1,
                        seed = 42,
                        font.label = c(12, "plain", "black"), label.rectangle = FALSE,
                        palette = c("#B31B21", "#1465AC", "darkgray"),
                        top = 15, select.top.method = c("padj", "fc"),
                        label.select = NULL,
                        main = NULL, xlab = "Log2 mean expression",  ylab = "Log2 fold change",
                        ggtheme = theme_classic(),...)
{
  
  if(!base::inherits(data, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if(!is.null(detection_call)){
    if(nrow(data)!=length(detection_call))
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if("detection_call" %in% colnames(data)){
    detection_call <- as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))
  
  # Legend position
  if(is.null(list(...)$legend)) legend <- c(0.12, 0.9)
  # If basemean logged, we'll leave it as is, otherwise log2 transform
  is.basemean.logged <- "baseMeanLog2" %in% colnames(data)
  if(is.basemean.logged){
    data$baseMean <- data$baseMeanLog2
  }
  else if("baseMean" %in% colnames(data)){
    data$baseMean <- log2(data$baseMean +1)
  }
  
  # Check data format
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(data))
  if(length(ss)>0) stop("The colnames of data must contain: ",
                        paste(ss, collapse = ", "))
  
  if(is.null(genenames)) genenames <- rownames(data)
  else if(length(genenames)!=nrow(data))
    stop("genenames should be of length nrow(data).")
  
  sig <- rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 1
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange,
                     padj = data$padj, sig = sig)
  
  # Change level labels
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(
    paste0("Up: ", sum(sig == 1)),
    paste0("Down: ", sum(sig == 2)),
    "NS"
  ) %>% .[.lev]
  
  data$sig <- factor(data$sig, labels = new.levels)
  
  
  # Ordering for selecting top gene
  select.top.method <- match.arg(select.top.method)
  if(select.top.method == "padj") data <- data[order(data$padj), ]
  else if(select.top.method == "fc") data <- data[order(abs(data$lfc), decreasing = TRUE), ]
  # select data for top genes
  complete_data <- stats::na.omit(data)
  labs_data <- subset(complete_data, padj <= fdr & name!="" & abs(lfc) >= log2(fc))
  labs_data <- utils::head(labs_data, top)
  # Select some specific labels to show
  if(!is.null(label.select)){
    selected_labels  <- complete_data %>%
      subset(complete_data$name  %in% label.select, drop = FALSE)
    labs_data <- dplyr::bind_rows(labs_data, selected_labels) %>%
      dplyr::distinct(.data$name, .keep_all = TRUE)
  }
  
  
  font.label <- list(size = as.numeric(font.label[1]), color=font.label[3], face=font.label[2])
  font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", font.label$face)
  
  # Plot
  NSandInSubfile = (data$name %in% rownames(subsetfile)) & (data$sig == "NS") #to change
  sigandinsubfile = (data$name %in% rownames(subsetfile)) & (data$sig != "NS")
  genes1 = subset(data, NSandInSubfile)
  genes2 = subset(data, sigandinsubfile)
  mean <- lfc <- sig <- name <- padj <-  NULL
  p <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig), size = size, alpha = alpha) +
    geom_point(data=genes1, color=subsetcolor1) + 
    geom_point(data=genes2, color=subsetcolor2, size=2.5) +
    geom_label_repel(data = genes2, aes(label = name), 
                     box.padding = unit(0.3, "lines"), 
                     segment.size = 0.5,
                     segment.color = subsetcolor2,
                     size = 3,
                     nudge_x = 1.5, nudge_y = 1, max.overlaps = Inf) + 
    geom_segment(data = genes2, aes(x = mean, y = lfc, xend = mean, yend = lfc),
                 arrow = arrow(length = unit(0.02, "npc")))
  
  max.overlaps = getOption("ggrepel.max.overlaps", default = Inf)
  
  if(label.rectangle){
    p <- p + ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = name),
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 1, seed = seed, fontface = font.label$face,
                                       size = font.label$size/3, color = font.label$color,
                                       max.overlaps = max.overlaps)
  }
  else{
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = name),
                                      box.padding = unit(0.35, "lines"),
                                      point.padding = unit(0.3, "lines"),
                                      force = 1, seed = seed, fontface = font.label$face,
                                      size = font.label$size/3, color = font.label$color,
                                      max.overlaps = max.overlaps)
  }
  
  p <- p + scale_x_continuous(breaks=seq(0, max(data$mean), 2))+
    labs(x = xlab, y = ylab, title = main, color = "")+ # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black"))
  
  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  p
}

#todo : add a helper function which is the git of ggmaplot, edit the sig color conditions and 
# add a condition that takes rownames(res) and color those that are FBF2 bound
# checking what is the limit
# change the proportions - zoom in
#make MA plots for both shrunken and regular
do_MA <-function(res, subfile, subcolor1, subcolor2) {
  myggmaplot(res, subfile, subcolor1, subcolor2, main = expression("MA plot"),
             fdr = 0.05, fc = 0, size = 1, palette = c("#B31B21", "#007FFF", "darkgray"),
             genenames = as.vector(res$name),
             ggtheme = ggplot2::theme_minimal(),
             top = 0)
}


cluster_expression <- function(dds, k, coldata, genes) {
  normalized_out <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5)) # get normalized count data from dds object
  mat<-assay(normalized_out)[rownames(genes), rownames(coldata)] # sig genes x samples 
  colnames(mat) <- rownames(coldata)
  mat.scaled <- t(apply(mat, 1, scale)) # center and scale each column (Z-score) then transpose
  colnames(mat.scaled)<-colnames(mat)
  
  h1 <- Heatmap(mat.scaled, cluster_rows = F, cluster_columns = F, km=k,
                column_labels = colnames(mat.scaled), name = "z_score", show_row_names = FALSE)
  
  print(h1)
}


cluster_mean_expression <- function(dds, k, coldata, genes, data_type = NULL, col = "slice") {
  # Perform variance stabilizing transformation
  normalized_out <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5))
  
  # Ensure that only genes present in the counts table are used
  valid_genes <- rownames(genes)[rownames(genes) %in% rownames(assay(normalized_out))]
  cat("Number of valid genes:", length(valid_genes), "\n")
  
  if (length(valid_genes) == 0) {
    cat("No valid genes found for the current data type. Exiting function.\n")
    return(NULL)
  }
  
  # Extract the relevant genes from the normalized matrix
  mat <- assay(normalized_out)[valid_genes, , drop = FALSE]
  
  # Initialize a list to store row means for each unique value in the column
  value_means <- list()
  
  # Loop through each unique value in the specified column
  for (val in unique(coldata[[col]])) {
    # Extract indices for the current value
    value_indices <- which(coldata[[col]] == val)
    
    # Combine worm, slice, and slice_index information to form the full column names
    full_column_names <- paste0(coldata$worm[value_indices], "_", coldata$slice[value_indices], "_", coldata$slice_index[value_indices])
    
    # Extract data for the current slices
    value_data <- mat[, full_column_names, drop = FALSE]
    
    # Calculate the row means for the current value
    value_mean <- rowMeans(value_data)
    
    # Store the row means for the current value
    col_name <- paste0("Value_", val)
    value_means[[col_name]] <- value_mean
  }
  
  # Combine row means for all values into a single matrix
  mat_means <- do.call(cbind, value_means)
  
  # Remove rows with zero variance and scale the matrix
  zero_var_rows <- apply(mat_means, 1, function(x) var(x) == 0)
  mat.scaled <- t(apply(mat_means[!zero_var_rows, ], 1, scale))
  
  # Set column names for the scaled matrix
  colnames(mat.scaled) <- names(value_means)
  
  # Generate heatmap with k-means clustering
  h1 <- Heatmap(mat.scaled, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                km = k, 
                use_raster = FALSE,  # Disable rasterization
                row_labels = NULL,  # Optional: hide row labels to fit more genes
                show_row_names = FALSE,  # Optional: hide row names for better clarity
                column_labels = colnames(mat.scaled), 
                name = "z_score", 
                height = unit(12, "cm"),  # Adjust height for better visualization
                column_title = paste("Heatmap for type:", data_type))
  
  h1_d <- draw(h1)
  return(h1_d)
}


plot_all_sRNA_type_heatmaps <- function(dds, coldata, gene_list, k = 3) {
  # Get the unique sRNA types from the gene list
  sRNA_types <- unique(gene_list$Type)
  
  # Loop through each sRNA type
  for (sRNA_type in sRNA_types) {
    cat("Processing sRNA type:", sRNA_type, "\n")
    
    # Filter gene list to get WB IDs for the current sRNA type
    genes_of_type <- gene_list$WB[gene_list$Type == sRNA_type]
    genes_of_type_df <- data.frame(row.names = genes_of_type, WB = genes_of_type)
    
    # Call the cluster_mean_expression function for the current sRNA type
    if (length(genes_of_type) > 0) {  # Ensure there are genes to process
      cluster_mean_expression(dds, k, coldata, genes_of_type_df, sRNA_type)
    } else {
      cat("No genes found for sRNA type:", sRNA_type, "\n")
    }
  }
}



plot_by_type <- function(cluster, annot, num, chr) {
  cluster_df <- data.frame(ID =  c(cluster))
  cluster_df$Type <- annot$Type[match(cluster_df$ID, annot$GeneID)]
  ggplot(cluster_df, aes(x = Type)) +
    geom_bar() +
    labs(x = "Type", y = "Count", title = glue('Count of Genes by Type - Cluster {num}, {chr}'))
}

generate_correlation_heatmap <- function(dds, coldata, col) {
  # Perform variance stabilizing transformation
  normalized_out <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5))
  
  # Initialize a list to store average counts for each unique value in the specified column
  value_means <- list()
  
  # Get unique values of the specified column
  unique_values <- unique(coldata[[col]])
  
  # Loop through each unique value
  for (val in unique_values) {
    # Extract indices for the current value
    value_indices <- which(coldata[[col]] == val)
    
    # Combine worm, slice, and slice_index information to form the new sample names
    full_column_names <- paste0(coldata$worm[value_indices], "_", coldata$slice[value_indices], "_", coldata$slice_index[value_indices])
    
    # Extract the corresponding columns from the normalized counts matrix
    value_data <- assay(normalized_out)[, full_column_names, drop = FALSE]
    
    # Calculate the row means across these columns for the current value
    value_mean <- rowMeans(value_data)
    
    # Store the row means for the current value
    value_means[[paste0("Value_", val)]] <- value_mean
  }
  
  # Combine row means for all values into a single matrix
  mat_means <- do.call(cbind, value_means)
  
  # Calculate Pearson correlation across values
  correlation_matrix <- cor(mat_means, method = "pearson")
  
  # Create heatmap using ComplexHeatmap
  heatmap_colors <- colorRamp2(c(0.7, 0.85, 1), c("blue", "white", "red"))
  Heatmap(correlation_matrix, 
          name = "Pearson Correlation", 
          col = heatmap_colors,
          cluster_rows = FALSE, 
          cluster_columns = FALSE,
          show_row_names = TRUE, 
          show_column_names = TRUE,
          column_labels = names(value_means),
          row_labels = names(value_means))
}



plot_param_correlation <- function(dds, coldata, col) {
  # Perform variance stabilizing transformation
  normalized_out <- vst(dds, blind = FALSE, nsub = sum(rowMeans(counts(dds, normalized = TRUE)) > 5))
  
  # Get unique values for the column specified
  unique_values <- unique(coldata[[col]])
  
  # Loop over each unique value
  for (val in unique_values) {
    # Extract indices for the current value
    indices <- which(coldata[[col]] == val)
    
    # Create updated sample names (worm_slice_slice_index) and match them with the colData
    full_column_names <- paste0(coldata$worm[indices], "_", coldata$slice[indices], "_", coldata$slice_index[indices])
    
    # Extract data for the current slices using the updated sample names
    data <- assay(normalized_out)[, full_column_names, drop = FALSE]
    
    # Calculate Pearson correlation for the extracted data
    correlation_matrix <- cor(data, method = "pearson")
    
    # Plot the correlation matrix as a heatmap
    heatmap_title <- paste("Correlation for ", val)
    heatmap <- Heatmap(correlation_matrix, 
                       name = "Correlation", 
                       col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
                       cluster_rows = FALSE, 
                       cluster_columns = FALSE,
                       show_row_names = TRUE, 
                       show_column_names = TRUE,
                       column_labels = colnames(data),
                       row_labels = colnames(data),
                       top_annotation = HeatmapAnnotation(slice = anno_text(coldata$slice[indices], rot = 0)),
                       column_title = heatmap_title)
    
    # Display the heatmap
    print(heatmap)
  }
}



calculate_missing_percentages <- function(gene_list, dds_sRNAs) {
  # Step 1: Filter genes on X and not on X for snRNA, snoRNA, and piRNA
  
  # On X chromosome
  X_genes_filtered <- gene_list %>%
    filter(Chr == "X" & Type %in% c("snRNA", "snoRNA", "piRNA"))
  
  # Not on X chromosome (autosomes)
  non_X_genes_filtered <- gene_list %>%
    filter(Chr != "X" & Type %in% c("snRNA", "snoRNA", "piRNA"))
  
  # Extract list of WB IDs from dds_sRNAs count table
  dds_genes <- rownames(counts(dds_sRNAs))
  
  # Initialize a list to store the results
  results <- list()
  
  # Step 2: Filter missing genes and calculate percentages and counts for each type
  
  for (gene_type in c("snRNA", "snoRNA", "piRNA")) {
    
    # Missing on X
    missing_on_X <- X_genes_filtered %>%
      filter(Type == gene_type & !(WB %in% dds_genes))
    total_on_X <- nrow(X_genes_filtered %>% filter(Type == gene_type))
    missing_on_X_count <- nrow(missing_on_X)
    missing_on_X_pct <- (missing_on_X_count / total_on_X) * 100
    
    # Missing not on X (autosomes)
    missing_non_X <- non_X_genes_filtered %>%
      filter(Type == gene_type & !(WB %in% dds_genes))
    total_non_X <- nrow(non_X_genes_filtered %>% filter(Type == gene_type))
    missing_non_X_count <- nrow(missing_non_X)
    missing_non_X_pct <- (missing_non_X_count / total_non_X) * 100
    
    # Store the results in the list
    results[[gene_type]] <- list(
      missing_on_X_count = missing_on_X_count,
      total_on_X = total_on_X,
      missing_on_X_pct = missing_on_X_pct,
      missing_non_X_count = missing_non_X_count,
      total_non_X = total_non_X,
      missing_non_X_pct = missing_non_X_pct
    )
  }
  
  # Step 3: Return the results
  return(results)
}

# Function to update sample names in count table and colData with slice_index
update_sample_names <- function(counts, coldata) {
  
  # Ensure that coldata has a column for slice_index and worm
  if (!("slice_index" %in% colnames(coldata)) || !("worm" %in% colnames(coldata)) || !("slice" %in% colnames(coldata))) {
    stop("colData must have 'slice_index', 'worm', and 'slice' columns.")
  }
  
  # Create new sample names by concatenating worm, slice, and slice_index
  coldata$new_sample_name <- paste0(coldata$worm, "_", coldata$slice, "_", coldata$slice_index)
  
  # Update the column names of the count table based on the new sample names
  colnames(counts) <- coldata$new_sample_name
  
  # Optionally update the rownames of colData to reflect the new sample names
  rownames(coldata) <- coldata$new_sample_name
  
  # Return updated counts and coldata
  return(list(updated_counts = counts, updated_coldata = coldata))
}


raj_heatmap <- function(count_data, col_data) {
  # Step 1: Order sample names based on percent_distal_to_proximal
  col_data <- col_data %>%
    arrange(percent_distal_to_proximal)
  ordered_samples <- rownames(col_data)
  ordered_samples <- intersect(ordered_samples, colnames(count_data))
  
  if (length(ordered_samples) == 0) {
    stop("No matching samples found between count data and design data.")
  }
  
  # Step 2: Subset and order count_data by ordered samples
  data_matrix <- as.matrix(count_data[, ordered_samples, drop = FALSE])
  
  # Step 3: Precompute row clustering based on Pearson correlation
  dist_matrix <- as.dist(1 - cor(t(data_matrix), method = "pearson"))
  row_cluster <- hclust(dist_matrix, method = "average")
  row_dend <- as.dendrogram(row_cluster)
  
  # Step 4: Use ComplexHeatmap’s color scaling to visualize deviations
  Heatmap(
    data_matrix,
    name = "Expression Level",
    cluster_rows = row_dend,           # Apply precomputed clustering
    cluster_columns = FALSE,            # Keep columns ordered as-is
    col = colorRamp2(c(-2, 0, 2), c("purple", "cyan", "yellow")),  # Z-score color range
    row_title = "Genes",
    column_title = "Samples",
    heatmap_legend_param = list(
      title = "Z-score",
      at = c(-2, 0, 2),
      labels = c("Below Mean -2SD", "Mean", "Above Mean +2SD")
    ),
    show_row_names = FALSE              # Adjust as needed for clarity
  )
}



#~~~~~~~~~~~~~~x_chr_proj main~~~~~~~~~~~~~~~~~~~
# gene_list <- read.delim("sRNAs.tsv", header = TRUE, sep = "\t")
# gene_list <- gene_list[!duplicated(gene_list$WB), ]
# # Extract gene IDs from the 'WB' column
# gene_ids <- gene_list$WB
# non_rRNA_genes <- gene_list$WB[gene_list$Type != "rRNA"]
# 
# colData <- prep_coldata("des_start_end.csv")
# counts <- prep_counts("counts_all.csv")
# 
# dds <- prepare_data(counts, "end", colData)
# do_pca(dds, 1, 2)
# 
# 
# dds_sRNAs <- dds[rownames(dds) %in% gene_ids, ]
# do_pca(dds_sRNAs, 1, 2)
# 
# 
# plot_slice_correlation(dds_sRNAs, colData)
# 
# outlier_samples <- c("W1_S2", "W2_S2", "W3_S3",
#                         "W1_S8", "W2_S8",  "W5_S13","W2_S14", "W1_S15", "W1_S16")
# 
# # Filter out the outlier samples from colData
# colData_filtered <- colData[!rownames(colData) %in% outlier_samples, ]
# colData_filtered <- colData_filtered[order(as.numeric(gsub("S", "", colData_filtered$slice))), ]
# 
# dds_sRNAs <- dds_sRNAs[, !colnames(dds_sRNAs) %in% outlier_samples]
# do_pca(dds_sRNAs, 1, 2)
# plot_slice_correlation(dds_sRNAs, colData_filtered)
# 
# do_pca_means(dds_sRNAs, colData_filtered)
# generate_correlation_heatmap(dds_sRNAs, colData_filtered)
# 
# plot_all_sRNA_type_heatmaps(dds_sRNAs, colData_filtered, gene_list, k = 4)
# 
# # comparing x to autosomes  
# X_genes <- gene_list$WB[gene_list$Chr == "X"]
# ddsonlyX <- dds_sRNAs[rownames(dds_sRNAs) %in% X_genes,]
# generate_correlation_heatmap(ddsonlyX, colData_filtered)
# 
# plot_all_sRNA_type_heatmaps(ddsonlyX, colData_filtered, gene_list, k = 3)
# 
# 
# ddsAuthosomes <- dds_sRNAs[!(rownames(dds_sRNAs) %in% X_genes),]
# generate_correlation_heatmap(ddsAuthosomes, colData_filtered)
# 
# plot_all_sRNA_type_heatmaps(ddsAuthosomes, colData_filtered, gene_list, k = 3)
# 
# # ~~~~~~~~~~~~~~~ checking missing genes for silencing ~~~~~~~~~~~~~~~ 
# 
# missing_percentages <- calculate_missing_percentages(gene_list, dds_sRNAs)
# 
# # Printing results
# for (gene_type in names(missing_percentages)) {
#   cat(gene_type, ":\n")
#   cat("  - Missing on X: ", missing_percentages[[gene_type]]$missing_on_X_count, 
#       " out of ", missing_percentages[[gene_type]]$total_on_X, 
#       " (", missing_percentages[[gene_type]]$missing_on_X_pct, "%)\n", sep="")
#   cat("  - Missing on autosomes: ", missing_percentages[[gene_type]]$missing_non_X_count, 
#       " out of ", missing_percentages[[gene_type]]$total_non_X, 
#       " (", missing_percentages[[gene_type]]$missing_non_X_pct, "%)\n", sep="")
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~ rajewsky miRNA table ~~~~~~~~~~~~~~~~~~~~~

colData_raj <- prep_coldata("design_miRNA_with_slices_dbscan.csv")
counts_raj <- prep_counts("Rajewsky_miRNA_count.csv")
cpm_raj <- prep_counts("miRNA_cpm_table.csv")

raj_heatmap(cpm_raj, colData_raj)

# Filter out samples with slice_index == -1 (outliers)
filtered_colData_raj = colData_raj[colData_raj['slice_index'] != -1,]

# Filter counts table to keep only columns that match the filtered colData
filtered_counts_raj <- counts_raj[, rownames(filtered_colData_raj)]

updated_data <- update_sample_names(filtered_counts_raj, filtered_colData_raj)

# Extract the updated counts and colData
updated_counts_raj <- updated_data$updated_counts
updated_colData_raj <- updated_data$updated_coldata

dds_raj <- prepare_data(updated_counts_raj, "end", updated_colData_raj)
do_pca(dds_raj, 1, 2)

plot_param_correlation(dds_raj, colData(dds_raj), "slice_index")

outlier_samples_raj <- c("W2_S13_13", "W5_S13_13", "W2_S14_14", "W1_S16_14", "W2_S2_2",  "W2_S3_3", "W5_S3_3", "W3_S5_5", "W3_S4_4", "W3_S6_6", "W3_S7_7", "W3_S8_8", "W3_S9_9")

# Filter out the outlier samples from colData_raj
colData_raj_filtered <- updated_colData_raj[!rownames(updated_colData_raj) %in% outlier_samples_raj, ]
colData_raj_filtered <- colData_raj_filtered[order(colData_raj_filtered$slice_index), ]
colData_raj_filtered$worm <- droplevels(colData_raj_filtered$worm)

# Filter out the outliers from dds_raj
dds_raj_filtered <- dds_raj[, !colnames(dds_raj) %in% outlier_samples_raj]
dds_raj_filtered <- dds_raj_filtered[, rownames(colData_raj_filtered)]
colData(dds_raj_filtered)$worm <- droplevels(colData(dds_raj_filtered)$worm)

plot_param_correlation(dds_raj_filtered, colData(dds_raj_filtered), "slice_index")
do_pca(dds_raj_filtered, 1, 2)

x_raj <- read.delim("mature_ids_chrX.tsv", header = TRUE, sep = "\t")
row.names(x_raj) <- x_raj$mature_id

ddsRajonlyX <- dds_raj_filtered[rownames(dds_raj_filtered) %in% x_raj$mature_id,]
do_pca(ddsRajonlyX, 1, 2)
generate_correlation_heatmap(ddsRajonlyX, colData(ddsRajonlyX),  "slice_index")
cluster_mean_expression(ddsRajonlyX, k = 3, coldata = colData(ddsRajonlyX), genes = x_raj, data_type = "X chr", col = "slice_index")
res_mir_X <-  results_data(ddsRajonlyX, "start", "end")
do_heatmap_means(res_mir_X, ddsRajonlyX, colData(ddsRajonlyX), c('start', 'end'), col = "slice_index")


ddsRajAutosomes <- dds_raj_filtered[!(rownames(dds_raj_filtered) %in% x_raj$mature_id), ]
# Create a data frame with rownames of ddsRajAutosomes as both rownames and values
df_autosomes <- data.frame(mature_id = rownames(ddsRajAutosomes), row.names = rownames(ddsRajAutosomes))
generate_correlation_heatmap(ddsRajAutosomes, colData(ddsRajAutosomes),  "slice_index")
cluster_mean_expression(ddsRajAutosomes, k = 3, coldata = colData(ddsRajAutosomes), genes = df_autosomes, data_type = "autosomes", col = "slice_index")
res_mir_auto <-  results_data(ddsRajAutosomes, "start", "end")
do_heatmap_means(res_mir_auto, ddsRajAutosomes, colData(ddsRajAutosomes), c('start', 'end'), col = "slice_index")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~ rajewsky mRNA table ~~~~~~~~~~~~~~~~~~~~~

colData_raj_mRNA <- prep_coldata("design_mRNA_with_slices_dbscan.csv")
counts_raj_mRNA <- prep_counts("raj_mRNA_table.csv")

# Filter out samples with slice_index == -1 (outliers)
filtered_colData_raj_mRNA = colData_raj_mRNA[colData_raj_mRNA['slice_index'] != -1,]

# Filter counts table to keep only columns that match the filtered colData
filtered_counts_raj_mRNA <- counts_raj_mRNA[, rownames(filtered_colData_raj_mRNA)]

updated_data_mRNA <- update_sample_names(filtered_counts_raj_mRNA, filtered_colData_raj_mRNA)

# Extract the updated counts and colData
updated_counts_raj_mRNA <- updated_data_mRNA$updated_counts
updated_colData_raj_mRNA <- updated_data_mRNA$updated_coldata

dds_raj_mRNA <- prepare_data(updated_counts_raj_mRNA, "end", updated_colData_raj_mRNA)

do_pca(dds_raj_mRNA, 1, 2)

plot_param_correlation(dds_raj_mRNA, colData(dds_raj_mRNA), "slice_index")

outlier_samples_raj <- c("W1_S1_1", "W4_S5_1", "W5_S1_1", "W6_S3_1", "W2_S4_2",
                         "W6_S10_8", "W5_S7_8", "W5_S8_9", "W6_S11_9", 
                         "W1_S2_2", "W4_S10_6", "W2_S10_9", "W1_S10_10", 
                         "W2_S11_10", "W2_S12_11", "W2_S13_12", "W1_S15_14")

# Filter out the outlier samples from colData_raj
colData_raj_mRNA_filtered <- updated_colData_raj_mRNA[!rownames(updated_colData_raj_mRNA) %in% outlier_samples_raj, ]
colData_raj_mRNA_filtered <- colData_raj_mRNA_filtered[order(colData_raj_mRNA_filtered$slice_index), ]
colData_raj_mRNA_filtered$worm <- droplevels(colData_raj_mRNA_filtered$worm)

# Filter out the outliers from dds_raj
dds_raj_mRNA_filtered <- dds_raj_mRNA[, !colnames(dds_raj_mRNA) %in% outlier_samples_raj]
dds_raj_mRNA_filtered <- dds_raj_mRNA_filtered[, rownames(colData_raj_mRNA_filtered)]
colData(dds_raj_mRNA_filtered)$worm <- droplevels(colData(dds_raj_mRNA_filtered)$worm)

plot_param_correlation(dds_raj_mRNA_filtered, colData(dds_raj_mRNA_filtered), "slice_index")
do_pca(dds_raj_mRNA_filtered, 1, 2)
generate_correlation_heatmap(dds_raj_mRNA_filtered, colData(dds_raj_mRNA_filtered),  "slice_index")

x_raj_mRNA <- read.delim("unique_chrX_transcripts.csv", header = TRUE, sep = "\t")
row.names(x_raj_mRNA) <- x_raj_mRNA$transcript_id

ddsRajmRNAonlyX <- dds_raj_mRNA_filtered[rownames(dds_raj_mRNA_filtered) %in% x_raj_mRNA$transcript_id,]
generate_correlation_heatmap(ddsRajmRNAonlyX, colData(ddsRajmRNAonlyX),  "slice_index")
cluster_mean_expression(ddsRajmRNAonlyX, k = 2, coldata = colData(ddsRajmRNAonlyX), genes = x_raj_mRNA, data_type = "X chr", col = "slice_index")
res_mRNA_X <-  results_data(ddsRajmRNAonlyX, "start", "end")
do_heatmap_means(res_mRNA_X, ddsRajmRNAonlyX, colData(ddsRajmRNAonlyX), c('start', 'end'), col = "slice_index")

ddsRajmRNAAutosomes <- dds_raj_mRNA_filtered[!(rownames(dds_raj_mRNA_filtered) %in% x_raj_mRNA$transcript_id), ]
# Create a data frame with rownames of ddsRajAutosomes as both rownames and values
df_autosomes <- data.frame(transcript_id = rownames(ddsRajmRNAAutosomes), row.names = rownames(ddsRajmRNAAutosomes))
generate_correlation_heatmap(ddsRajmRNAAutosomes, colData(ddsRajmRNAAutosomes),  "slice_index")
cluster_mean_expression(ddsRajmRNAAutosomes, k = 2, coldata = colData(ddsRajmRNAAutosomes), genes = df_autosomes, data_type = "autosomes", col = "slice_index")
res_mRNA_auto <-  results_data(ddsRajmRNAAutosomes, "start", "end")
do_heatmap_means(res_mRNA_auto, ddsRajmRNAAutosomes, colData(ddsRajmRNAAutosomes), c('start', 'end'), col = "slice_index")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ordered_coldata <- colData[colnames(dds), ]
res_start <-  results_data(dds, "start", "end")
res_intermediate <- results_data(dds, "intermediate", "end")

res_X_start_end <- res_start[rownames(res_start) %in% X_genes$genes,]
print(paste("up between start end on X:", count_upregulated(res_X_start_end)))
print(paste("down between start end on X:", count_downregulated(res_X_start_end)))
do_heatmap(res_X_start_end, ddsonlyX, ordered_coldata, c('start', 'end'))
do_heatmap_means(res_X_start_end, ddsonlyX, ordered_coldata, c('start', 'end'))

up_x <- extract_upregulated_genes(res_X_start_end, "start", "end")
down_x <- extract_downregulated_genes(res_X_start_end, "start", "end")

res_autho_start_end <- res_start[!(rownames(res_start) %in% X_genes$genes),]
print(paste("up between start end on autosomes:", count_upregulated(res_autho_start_end)))
print(paste("down between start end on autosomes:", count_downregulated(res_autho_start_end)))
do_heatmap(res_autho_start_end, ddsAuthosomes, ordered_coldata, c('start', 'end'))
do_heatmap_means(res_autho_start_end, ddsAuthosomes, ordered_coldata, c('start', 'end'))


up_auto <- extract_upregulated_genes(res_autho_start_end, "start", "end")
down_auto <- extract_downregulated_genes(res_autho_start_end, "start", "end")

res_X_inter_end <- res_intermediate[rownames(res_intermediate) %in% X_genes$genes,]
print(paste("up between intermediate end on X:", count_upregulated(res_X_inter_end)))
print(paste("down between intermediate end on X:", count_downregulated(res_X_inter_end)))
do_heatmap(res_X_inter_end, ddsonlyX, ordered_coldata, c('intermediate', 'end'))

res_autho_inter_end <- res_intermediate[!(rownames(res_intermediate) %in% X_genes$genes),]
print(paste("up between intermediate end on autosomes:", count_upregulated(res_autho_inter_end)))
print(paste("down between intermediate end on autosomes:", count_downregulated(res_autho_inter_end)))
do_heatmap(res_autho_inter_end, ddsAuthosomes, ordered_coldata, c('intermediate', 'end'))


# stable genes
stable_start <- results_stable_genes(dds, "start", "end", 1.7)
stable <- stable_start[(stable_start$padj < 0.05) & (stable_start$baseMean) > 10,]
do_MA(res_start, stable, "#1dce00", "#ec7e0a")

stable_start_X <- stable[rownames(stable) %in% X_genes$genes,]
do_MA(res_X_start_end, stable_start_X, "#1dce00", "#ec7e0a")


stable_start_auto <- stable[!(rownames(stable) %in% X_genes$genes),]
do_MA(res_autho_start_end, stable_start_auto, "#1dce00", "#ec7e0a")

stable_up_down_x <- rbind(stable_start_X, up_x, down_x)
cluster_expression(dds, 5, ordered_coldata, stable_up_down_x)
h_x <- cluster_mean_expression(dds, 4, ordered_coldata, stable_up_down_x)

# extract clusters on x
x_clusters <- heatmap_extract_cluster(h_x, stable_up_down_x, which = "row")
cluster_1_x <- x_clusters[x_clusters$Cluster == 1,]$ID
cluster_2_x <- x_clusters[x_clusters$Cluster == 2,]$ID
cluster_3_x <- x_clusters[x_clusters$Cluster == 3,]$ID
cluster_4_x <- x_clusters[x_clusters$Cluster == 4,]$ID

annot <- read.table("sRNA_celegans_annotations.tsv", header = TRUE, sep = "\t")
annot_x <- annot[annot$Chr == "X",]

x_term <- annot_x %>% dplyr::distinct(Type, GeneID) %>% as.data.frame()

cluster1_enriched <- enricher(gene = cluster_1_x, TERM2GENE = x_term)
as.data.frame(cluster1_enriched)
plot_by_type(cluster_1_x, annot_x, 1, "X")

cluster2_enriched <- enricher(gene = cluster_2_x, TERM2GENE = x_term)
as.data.frame(cluster2_enriched)
fit_cluster2_x <- plot(barplot(cluster2_enriched, showCategory = 10))
plot_by_type(cluster_2_x, annot_x, 2, "X")

cluster3_enriched <- enricher(gene = cluster_3_x, TERM2GENE = x_term)
as.data.frame(cluster3_enriched)
fit_cluster3_x <- plot(barplot(cluster3_enriched, showCategory = 10))
plot_by_type(cluster_3_x, annot_x, 3, "X")

cluster4_enriched <- enricher(gene = cluster_4_x, TERM2GENE = x_term)
as.data.frame(cluster4_enriched)
fit_cluster4_x <- plot(barplot(cluster4_enriched, showCategory = 10))
plot_by_type(cluster_4_x, annot_x, 4, "X")


# extract clusters on autosomes
stable_up_down_auto <- rbind(stable_start_auto, up_auto, down_auto)
cluster_expression(dds, 3, ordered_coldata, stable_up_down_auto)
h_auto <- cluster_mean_expression(dds, 3, ordered_coldata, stable_up_down_auto)

#extract clusters on autosomes
annot_auto <- annot[!(annot$Chr == "X"),]
auto_term <- annot_auto %>% dplyr::distinct(Type, GeneID) %>% as.data.frame()

auto_clusters <- heatmap_extract_cluster(h_auto, stable_up_down_auto, which = "row")
cluster_1_auto <- auto_clusters[auto_clusters$Cluster == 1,]$ID
cluster_2_auto <- auto_clusters[auto_clusters$Cluster == 2,]$ID
cluster_3_auto <- auto_clusters[auto_clusters$Cluster == 3,]$ID

cluster1_enriched_auto <- enricher(gene = cluster_1_auto, TERM2GENE = auto_term)
as.data.frame(cluster1_enriched_auto)
fit_cluster1_auto <- plot(barplot(cluster1_enriched_auto, showCategory = 10))
plot_by_type(cluster_1_auto, auto_term, 1, "autosomes")

cluster2_enriched_auto <- enricher(gene = cluster_2_auto, TERM2GENE = auto_term)
as.data.frame(cluster2_enriched_auto)
fit_cluster2_auto <- plot(barplot(cluster2_enriched_auto, showCategory = 10))
plot_by_type(cluster_2_auto, auto_term, 2, "autosomes")

cluster3_enriched_auto <- enricher(gene = cluster_3_auto, TERM2GENE = auto_term)
as.data.frame(cluster3_enriched_auto)
fit_cluster3_auto <- plot(barplot(cluster3_enriched_auto, showCategory = 10))
plot_by_type(cluster_3_auto, auto_term, 3, "autosomes")
