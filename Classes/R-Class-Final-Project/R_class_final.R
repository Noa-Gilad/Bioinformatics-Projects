# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Enviorment Settings~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load necessary libraries for data analysis, visualization, and gene enrichment.
library(DESeq2)            # For differential expression analysis
library(ggplot2)           # For plots
library(ComplexHeatmap)    # For generating heatmaps
library(GeneOverlap)       # For gene overlap analysis
library(factoextra)        # For PCA visualization
library(ggforce)           # For enhanced ggplot features
library(ggpubr)            # For combining multiple plots
library(tibble)            # For managing data frames
library(clusterProfiler)   # For gene ontology enrichment analysis
library(org.Ce.eg.db)      # WormBase annotations for C. elegans
library(AnnotationDbi)     # Database interface for annotations
library(RColorBrewer)      # Color palettes for plots
library(circlize)          # Heatmap visualization utilities
library(extrafont)         # Extra fonts for plotting
library(rJava)             # Java support (required by venneuler)
library(venneuler)         # For generating proportional Venn diagrams
data(GeneOverlap)          # Load GeneOverlap dataset

# Set working directory to access necessary files.
setwd("C:\\Users\\gilad\\OneDrive\\Documents")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Constants~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set thresholds for filtering data based on read counts and statistical significance.
readsLim <- 10             # Minimum number of reads required for a gene to be included
lockBinding("readsLim", globalenv())  # Prevent modification of readsLim

alphaLim <- 0.05           # Statistical significance threshold for adjusted p-values
lockBinding("alphaLim", globalenv())  # Prevent modification of alphaLim

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepares the count table by filtering out low-count genes and removing specific unwanted rows.
prep_counts <- function(file) {
  counts <- read.delim(file, header=TRUE, row.names = 1, sep = ",")  # Load the count data from file
  row_names_to_remove <- c("WBGene00219609", "WBGene00219306", "WBGene00189951")  # Unwanted genes to remove
  counts <- counts[!(row.names(counts) %in% row_names_to_remove),]  # Filter out unwanted genes
  return(counts[which(rowSums(counts) >= readsLim),])  # Keep genes with total counts above threshold
}

# Prepare condition factors for DESeq2 analysis (experimental conditions: "Triple" and "WT").
prep_condition <- function() {
  return(factor(c("Triple", "Triple", "Triple", "WT", "WT", "WT")))  # Define conditions for each sample
}

# Prepare the column data (sample information) for DESeq2 based on counts and conditions.
prep_coldata <- function(counts, condition) {
  return(data.frame(row.names = colnames(counts), condition))  # Prepare column data for DESeq2 analysis
}

# Load the data, create a DESeq2 dataset, and perform differential expression analysis.
prepare_data <- function(counts, refer) {
  condition <- prep_condition()  # Get the condition (experimental grouping)
  coldata <- prep_coldata(counts, condition)  # Prepare column data for DESeq2
  
  # Construct DESeq2 dataset from count data and coldata
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
  dds$condition <- relevel(dds$condition, ref = refer)  # Set reference level for condition factor
  dds <- DESeq(dds)  # Perform DESeq2 analysis
  return(dds)  # Return DESeq2 dataset
}

# Perform PCA (Principal Component Analysis) on rlog-transformed data. 
# A modified versions of the opensource code of built-in deseq2 functions for better visualization
do_pca <- function(dds, pc1, pc2, intgroup="condition", ntop=500, returnData=FALSE) {
  transformed <- rlog(dds, blind=FALSE)  # Perform rlog transformation to stabilize variance
  
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
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC",pc1," (",round(percentVar[pc1] * 100),"%)")) +
    ylab(paste0("PC",pc2," (",round(percentVar[pc2] * 100),"%)")) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.title = element_blank(), text=element_text(size=12,  family="sans")) + 
    coord_fixed(xlim=c(-15, 15), ylim = c(-15, 15)) +
    ggforce::geom_mark_ellipse(aes(color = group, fill=group)) +  # Draw ellipses around groups
    scale_color_manual(values = c("blue","black")) +  # Customize colors
    scale_fill_manual(values = c("blue","black"))
}

# Extracts and sorts DESeq2 results, returning a summary and top results.
results_data <- function(dds, toCompare, baseLevel) {
  # Get DESeq2 results based on the comparison between 'toCompare' and 'baseLevel'
  res <- results(dds, alpha = alphaLim, contrast = c("condition", toCompare, baseLevel))
  
  # Display the top rows of the results in tidy format (for checking purposes)
  head(results(dds, tidy=TRUE))
  
  # Sort results by adjusted p-value
  res <- res[order(res$padj),]
  
  # Display the top 10 results (optional)
  head(res, 10)
  
  # Provide a summary of the DESeq2 results (e.g., number of significant genes)
  summary(res)
  
  # Return the sorted DESeq2 results
  return(res)
}

# Generates an MA plot with specific coloring for different gene sets
# A modified versions of the opensource code of ggmaplot function from ggpubr package for better visualization
do_MA <- function(data, subsetfile, subsetcolor1, subsetcolor2, subName = "FBF2", fdr = 0.05, fc = 0, genenames = as.vector(data$name),
                  detection_call = NULL, size = 1, alpha = 1, seed = 42,
                  font.label = c(12, "plain", "black"), label.rectangle = FALSE,
                  palette = c("red3", "dodgerblue3", subsetcolor1, subsetcolor2, "darkgray"),
                  top = 0, select.top.method = c("padj", "fc"), label.select = NULL,
                  main = expression("MA plot"), xlab = "Log2 mean expression", ylab = "Log2 fold change",
                  ggtheme = ggplot2::theme_minimal(), ...) {
  
  # Ensure the input data is of a valid class (matrix, data frame, or DESeq2 results)
  if(!base::inherits(data, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  
  # Check and handle detection calls (optional, defaults to '1' for all rows)
  if(!is.null(detection_call)){
    if(nrow(data) != length(detection_call))
      stop("detection_call must be a numeric vector of length = nrow(data)")
  } else if("detection_call" %in% colnames(data)){
    detection_call <- as.vector(data$detection_call)
  } else detection_call = rep(1, nrow(data))  # Default detection call is 1 for all rows
  
  # Set default legend position if not specified in additional arguments
  if(is.null(list(...)$legend)) legend <- c(0.12, 0.9)
  
  # Check if baseMean is already logged (log2 transformed), if not, log2 transform the baseMean column
  is.basemean.logged <- "baseMeanLog2" %in% colnames(data)
  if(is.basemean.logged){
    data$baseMean <- data$baseMeanLog2
  } else if("baseMean" %in% colnames(data)){
    data$baseMean <- log2(data$baseMean + 1)
  }
  
  # Ensure the required columns (baseMean, log2FoldChange, padj) are present in the data
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(data))
  if(length(ss) > 0) stop("The colnames of data must contain: ", paste(ss, collapse = ", "))
  
  # Handle gene names; if not provided, use row names
  if(is.null(genenames)) genenames <- rownames(data)
  else if(length(genenames) != nrow(data))
    stop("genenames should be of length nrow(data).")
  
  # Initialize a vector 'sig' to classify genes as upregulated, downregulated, or non-significant
  sig <- rep(3, nrow(data))
  
  # Classify downregulated genes based on FDR and fold-change thresholds
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call == 1)] = 2
  down <- sum(sig == 2)  # Count number of downregulated genes
  
  # Classify upregulated genes
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call == 1)] = 1
  up <- sum(sig == 1)  # Count number of upregulated genes
  
  # Label FBF-bound genes and include them in the plot's legend
  sig[rownames(data) %in% subsetfile$Suggested.Match & sig == 3] = 4  # Non-significant FBF-bound genes
  sig[rownames(data) %in% subsetfile$Suggested.Match & sig %in% c(1, 2)] = 5  # Significant FBF-bound genes
  
  # Create a data frame for plotting
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange, padj = data$padj, sig = sig)
  
  # Define factor levels for plot coloring (1: Upregulated, 2: Downregulated, 3: Non-significant, 4/5: FBF-bound)
  data$sig <- factor(data$sig, levels = c(1, 2, 4, 5, 3))
  
  # Customize the legend labels for the different significance levels
  new.levels <- c(
    paste0("Up: ", up),
    paste0("Down: ", down),
    paste0(subName, " and NS: ", sum(sig == 4)),  # Non-significant FBF-bound genes
    paste0(subName, " and significant: ", sum(sig == 5)),  # Significant FBF-bound genes
    "NS"  # Non-significant genes
  )
  
  # Update the factor levels with new legend labels
  former_data <- data
  data$sig <- factor(data$sig, labels = new.levels, levels = c(1, 2, 4, 5, 3))
  
  # Sort data based on padj or fold change (depending on 'select.top.method')
  select.top.method <- match.arg(select.top.method)
  if(select.top.method == "padj") data <- data[order(data$padj), ]
  else if(select.top.method == "fc") data <- data[order(abs(data$lfc), decreasing = TRUE), ]
  
  # Select top genes for labeling based on fold-change and p-value cutoffs
  complete_data <- stats::na.omit(data)
  labs_data <- subset(complete_data, padj <= fdr & name != "" & abs(lfc) >= log2(fc))
  labs_data <- utils::head(labs_data, top)
  
  # If specific genes are provided for labeling, include them in the labels
  if(!is.null(label.select)){
    selected_labels <- complete_data %>%
      subset(complete_data$name %in% label.select, drop = FALSE)
    labs_data <- dplyr::bind_rows(labs_data, selected_labels) %>%
      dplyr::distinct(.data$name, .keep_all = TRUE)
  }
  
  # Customize font settings for gene labels in the plot
  font.label <- list(size = as.numeric(font.label[1]), color=font.label[3], face=font.label[2])
  font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", font.label$face)
  
  # Subset FBF-bound genes for plotting
  NSandFBFbound <- (former_data$sig == 4)
  sigandFBFbound <- (former_data$sig == 5)
  genes1 <- subset(data, NSandFBFbound)  # Non-significant FBF-bound genes
  genes2 <- subset(data, sigandFBFbound)  # Significant FBF-bound genes
  
  # Create an MA plot using ggplot2
  p <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig), size = size, alpha = alpha) +  # Plot all points
    geom_point(data = genes1, color = subsetcolor1) +  # Highlight non-significant FBF-bound genes
    geom_point(data = genes2, color = subsetcolor2)  # Highlight significant FBF-bound genes
  
  # Customize label placement on the plot using ggrepel for avoiding label overlap
  max.overlaps <- getOption("ggrepel.max.overlaps", default = Inf)
  
  if(label.rectangle){
    # Label with rectangles around gene names
    p <- p + ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = name),
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 1, seed = seed, fontface = font.label$face,
                                       size = font.label$size/3, color = font.label$color,
                                       max.overlaps = max.overlaps)
  } else {
    # Label without rectangles around gene names
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = name),
                                      box.padding = unit(0.35, "lines"),
                                      point.padding = unit(0.3, "lines"),
                                      force = 1, seed = seed, fontface = font.label$face,
                                      size = font.label$size/3, color = font.label$color,
                                      max.overlaps = max.overlaps)
  }
  
  # Add x-axis breaks, labels, title, and other customizations to the plot
  p <- p + scale_x_continuous(breaks=seq(0, max(data$mean), 2)) +
    labs(x = xlab, y = ylab, title = main, color = "") +  # Title and axis labels
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black"))  # Horizontal lines at 0, -fc, and +fc
  
  # Apply ggplot2 theme and color palette
  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  
  # Return the final plot
  return(p)
}

# Creates a heatmap with customized annotations for FBF and germline enrichment.
do_heatmap <- function(res, dds, coldata, fileFBF, fileTitle, fileGermline) {
  sigs <- na.omit(res)  # Remove any NA values from the results
  sigs <- sigs[sigs$padj < alphaLim, ]  # Filter genes based on the adjusted p-value threshold (alphaLim)
  sigs.df <- as.data.frame(sigs)  # Convert results to a data frame
  sigs.df <- sigs.df[(abs(sigs.df$log2FoldChange) > 0), ]  # Keep only genes with significant log2 fold changes
  sigs.df <- sigs.df[order(sigs.df$log2FoldChange, decreasing = TRUE), ]  # Sort genes by log2 fold change in decreasing order
  
  # Add an indicator column for germline-enriched genes (1 if present in fileGermline, 0 otherwise)
  sigs.df$isGermlIne[rownames(sigs.df) %in% fileGermline$Suggested.Match] <- 1
  sigs.df$isGermlIne[!(rownames(sigs.df) %in% fileGermline$Suggested.Match)] <- 0
  
  # Add an indicator column for FBF-bound genes (1 if present in fileFBF, 0 otherwise)
  sigs.df$isFBFBound[rownames(sigs.df) %in% fileFBF$Suggested.Match] <- 1
  sigs.df$isFBFBound[!(rownames(sigs.df) %in% fileFBF$Suggested.Match)] <- 0
  
  # Perform regularized log transformation on the DESeq2 dataset to stabilize variance across samples
  rlog_out <- rlog(dds, blind=FALSE)
  
  # Extract the rlog-transformed expression matrix for the significant genes
  mat <- assay(rlog_out)[rownames(sigs.df), rownames(coldata)]
  colnames(mat) <- rownames(coldata)  # Assign sample names as column names
  
  # Scale the expression matrix by row (z-score normalization)
  mat.scaled <- t(apply(mat, 1, scale))  # Scale each row (gene)
  colnames(mat.scaled) <- colnames(mat)  # Restore column names after scaling
  
  # Convert germline-enriched gene indicators into a matrix
  isGermlIne <- as.matrix(sigs.df$isGermlIne)
  colnames(isGermlIne) <- "Is germlIne"  # Label the germline annotation
  
  # Convert FBF-bound gene indicators into a matrix
  isFBFBound <- as.matrix(sigs.df$isFBFBound)
  colnames(isFBFBound) <- fileTitle  # Label the FBF-bound annotation
  
  # Define the color scheme for the heatmap (blue for low expression, white for medium, red for high)
  colHeat <- colorRamp2(c(-1, 0, 1), c("blue4", "white", "red3"))
  
  # Create the main heatmap of scaled gene expression values
  h1 <- Heatmap(mat.scaled, cluster_rows = T, cluster_columns = T, 
                column_labels = colnames(mat.scaled), name = "z_score", col = colHeat)
  
  # for germline heatmap - put h2 in comment
  # Create a heatmap for germline enrichment annotation (binary: 1 for enriched, 0 for not)
  h2 <- Heatmap(isGermlIne, cluster_rows = F, name="is Germline enriched", col=c("khaki", "forestgreen"))
  
  # Create a heatmap for FBF-bound annotation (binary: 1 for FBF-bound, 0 for not)
  h3 <- Heatmap(isFBFBound, cluster_rows = F, name = fileTitle, col = c("gray", "maroon"))
  
  # Combine the heatmaps: gene expression, germline enrichment, and FBF-bound annotations
  h <- h1 + h2 + h3
  
  # for germline only heatmap do h<-h1+h3 and put the former in comment
  # h<-h1+h3
  
  # Display the combined heatmap
  print(h)
}

# Counts differentialy expressed genes (up or down) that are BFB bound
count_FBF_genes <- function(res, FBfile, direction = "upregulated") {
  # Extract upregulated or downregulated genes based on the direction
  genes <- extract_genes(res, ifelse(direction == "up", "up", "down"))
  
  # Count how many of the extracted genes are in the FBF-bound gene list
  FBF_genes_count <- sum(rownames(genes) %in% FBfile$Suggested.Match)
  
  return(FBF_genes_count)
}


# Extracts differentialy expressed genes (up or down) that are BFB bound
extract_FBF_genes <- function(res, FBfile, direction = "up") {
  # Extract upregulated or downregulated genes
  genes <- extract_genes(res, direction)
  
  # Filter the genes for FBF-bound genes
  FBF_genes <- genes[rownames(genes) %in% FBfile$Suggested.Match,]
  
  return(rownames(FBF_genes))
}

# Extracts differentialy expressed genes (up or down)
extract_genes <- function(res, direction = "up") {
  res <- na.omit(res)  # Remove rows with missing values (NA)
  
  # Filter for upregulated or downregulated genes based on the direction
  if (direction == "up") {
    genes <- res[res$log2FoldChange > 0 & res$padj < 0.05, ]
  } else {
    genes <- res[res$log2FoldChange < 0 & res$padj < 0.05, ]
  }
  
  # Add a new column for the gene names (WormBase IDs)
  genes$WB <- rownames(genes)
  
  return(genes)  # Return the rownames (WormBase IDs) of the filtered genes
}

# Prepares gene lists from a file for further analysis.
prepare_gene_lists <- function(file) {
  list <- read.delim(file, header = TRUE, row.names = NULL, sep = ",")  # Read the gene list from a file
  
  # Ensure column names are properly assigned
  colnames(list) <- colnames(list)[1:ncol(list)]
  
  # Remove entries with "not found" in the Suggested.Match column
  list <- list[!grepl("not found", list$Suggested.Match),]
  
  # Remove rows where Suggested.Match is NA or empty
  list <- list[!(is.na(list$Suggested.Match) | list$Suggested.Match == ""),]
  
  return(list)  # Return the cleaned gene list
}

# Prepares FBF1-bound gene lists, integrating with another file for WormBase IDs.
prepare_FBF1_list <- function(file, file2) {
  # Read the FBF1-bound gene list from file
  list1 <- read.delim(file, header = TRUE, row.names = NULL, sep = ",")
  
  # Read the gene name conversion file (Gene.name to Suggested.Match)
  list2 <- read.delim(file2, header = TRUE, row.names = NULL, sep = ",")
  
  # Match the Suggested.Match field in list1 with list2 based on Gene.name
  list1$Suggested.Match <- list2$Suggested.Match[match(list1$Gene.name, list2$Gene.name)]
  
  # Ensure column names are properly assigned
  colnames(list1) <- colnames(list1)[1:ncol(list1)]
  
  # Remove entries with "not found" in the Suggested.Match column
  list1 <- list1[!grepl("not found", list1$Suggested.Match),]
  
  # Remove rows where Suggested.Match is NA or empty
  list1 <- list1[!(is.na(list1$Suggested.Match) | list1$Suggested.Match == ""),]
  
  return(list1)  # Return the cleaned FBF1-bound gene list
}

# Generates a proportional Venn diagram with numbers
plot_overlap <- function(overlapObj, title) {
  
  # Extract the contingency table from the overlap test.
  # The table shows the number of genes in each category (in_A_not_in_B, not_in_A_in_B, etc.).
  contingency_table <- getContbl(overlapObj)
  
  # Assign the values from the contingency table to meaningful variables.
  not_in_A_not_in_B <- contingency_table[1, 1]
  in_A_not_in_B <- contingency_table[1, 2]
  not_in_A_in_B <- contingency_table[2, 1]
  in_A_and_B <- contingency_table[2, 2]  # This is the actual overlap between A and B
  
  # Extract additional statistics from the overlap object.
  pvalue <- getPval(overlapObj)  # P-value of the overlap
  odds_ratio <- getOddsRatio(overlapObj)  # Odds ratio of the overlap
  jaccard_index <- getJaccard(overlapObj)  # Jaccard index (measure of similarity between sets)
  
  # Create a proportional Venn diagram using the venneuler package.
  # The sizes of the circles (A and B) correspond to the number of genes in each set.
  venn_data <- venneuler(c(
    A = in_A_not_in_B + in_A_and_B,  # Genes in A (regulated genes)
    B = not_in_A_in_B + in_A_and_B,  # Genes in B (FBF-bound genes)
    "A&B" = in_A_and_B  # Overlap between A and B
  ))
  
  # Plot the Venn diagram with the title.
  plot(venn_data, main = title)
  
  # Annotate the numbers directly on the diagram (e.g., how many genes are in A, B, and the overlap).
  centers <- venn_data$centers  # Get the centers of the Venn diagram circles
  
  if (nrow(centers) >= 1) {
    # Annotate the number of genes in set A (regulated genes)
    grid.text(
      paste0("A: ", in_A_not_in_B + in_A_and_B), 
      x = centers[1, "x"] - 0.02, 
      y = centers[1, "y"], 
      gp = gpar(fontsize = 14, col = "black")
    )
  }
  
  if (nrow(centers) >= 2) {
    # Annotate the number of genes in set B (FBF-bound genes)
    grid.text(
      paste0("B: ", not_in_A_in_B + in_A_and_B), 
      x = centers[2, "x"] + 0.08, 
      y = centers[2, "y"], 
      gp = gpar(fontsize = 14, col = "black")
    )
  }
  
  # Annotate the overlap (number of genes in both A and B)
  if (in_A_and_B > 0) {
    grid.text(
      paste0("Overlap: ", in_A_and_B), 
      x = mean(centers[, "x"] + 0.06),  # Place annotation at the center of the overlap
      y = mean(centers[, "y"]), 
      gp = gpar(fontsize = 14, col = "black")
    )
  }
  
  # Annotate additional details beneath the diagram (p-value, odds ratio, Jaccard index, etc.)
  grid.text(
    paste("Not in A and B: ", not_in_A_not_in_B,  # Number of genes not in A and not in B
          "\nP-value =", format(pvalue, digits = 3),  # Format the p-value
          "\nOdds Ratio =", format(odds_ratio, digits = 3),  # Format the odds ratio
          "\nJaccard Index =", format(jaccard_index, digits = 3)),  # Format the Jaccard index
    x = 0.5, y = 0.1, gp = gpar(fontsize = 12)  # Place the annotation at the bottom of the plot
  )
  
}

# Performs the hypergeometric test and presents the proper venn diagram
test_overlap <- function(geneVecRegulated, FbfList, dds, refer, comp, fbf, regulation) {
  
  # Create the GeneOverlap object
  # This object compares two sets of genes (geneVecRegulated and FBF-bound genes from FbfList).
  # 'genome.size' is set to the total number of genes in the DESeq2 dataset (dds).
  overlapObj <- newGeneOverlap(geneVecRegulated$WB,
                               FbfList$Suggested.Match,
                               genome.size = length(dds))
  
  # Perform the overlap test (hypergeometric test or Fisher's exact test).
  # This determines whether the overlap between geneVecRegulated and FbfList is significant.
  overlapObj <- testGeneOverlap(overlapObj)
  print(overlapObj)  # Output the test results
  
  # Create a title for the Venn diagram based on the input parameters.
  # This title indicates the comparison being made (refer vs. comp), along with regulation and FBF.
  title <- paste(refer, "vs", comp, "- Fisher's exact test, A =", regulation, ", B =", fbf)
  
  plot_overlap(overlapObj, title)
}

# GO enrichment function with bar plot
perform_GO_enrichment <- function(wb_ids, ont = "BP", show_category = 15) {
  
  # Map WormBase gene IDs (wb_ids) to gene symbols using the org.Ce.eg.db database
  # "SYMBOL" refers to the common gene name, "WORMBASE" refers to the input key type (WormBase IDs)
  gene_symbols <- mapIds(org.Ce.eg.db, 
                         keys = wb_ids,  # WormBase gene IDs as input keys
                         column = "SYMBOL",  # Output gene symbols (common names)
                         keytype = "WORMBASE",  # Input key type is WormBase IDs
                         multiVals = "first")  # If multiple mappings, take the first one
  
  # Remove any NA values from the gene symbols (genes without a corresponding symbol)
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
  # Perform GO enrichment analysis using clusterProfiler's enrichGO function
  go_results <- enrichGO(
    gene         = gene_symbols,  # List of gene symbols to analyze
    OrgDb        = org.Ce.eg.db,  # C. elegans organism database for GO terms
    keyType      = "SYMBOL",  # Input gene IDs are gene symbols
    ont          = ont,  # Ontology category: "BP" (Biological Process) by default
    pAdjustMethod = "BH",  # Method for adjusting p-values (Benjamini-Hochberg)
    pvalueCutoff = 0.05,  # P-value cutoff for significant terms
    qvalueCutoff = 0.2  # q-value cutoff (FDR-adjusted p-value) for filtering GO terms
  )
  
  # Convert the GO enrichment results to a data frame for easier manipulation and viewing
  go_results_df <- as.data.frame(go_results)
  
  # Create a bar plot for the top GO categories (default is 15 categories)
  go_plot <- barplot(go_results, showCategory = show_category)  # Number of categories to show in the plot
  
  # Return a list containing the results data frame and the bar plot
  return(list(
    results_df = go_results_df,  # Data frame of GO enrichment results
    plot = go_plot  # Bar plot of the top GO categories
  ))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~main~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare lists of known FBF2-bound, FBF1-bound, and germline-enriched genes
FBF2bound <- prepare_gene_lists("FBF-2_pull_down_RNA.csv")
FBF1bound <- prepare_FBF1_list("FBF-1_pull_down_RNA.csv", "FBF-1_conversion.csv")
germline <- read.delim("list_of_germline_enriched_genes_WB_ID _(reinke_list).csv", header=TRUE, row.names = NULL, sep = ",")
colnames(germline) <- c("Suggested.Match")

#~~~~~~~~~~~~~~WT_reference~~~~~~~~~~~~~~~~~~~
# Prepare count data and create DESeq2 dataset
coldata1 <- prep_coldata(prep_counts("WT_vs_Triple_analysis.csv"), prep_condition())
counts <- prep_counts("WT_vs_Triple_analysis.csv")
dds1 <- prepare_data(counts, "WT")

# Perform PCA analysis
do_pca(dds1, 1, 2)

#~~~~~~~~~~~~~~~~~WT_VS_Triple~~~~~~~~~~~~~~~~~~~
# Extract DESeq2 results for WT vs Triple comparison
res_wt_triple <- results_data(dds1, "Triple", "WT")

# Filter results for germline-enriched genes
germlineEnriched <- (rownames(res_wt_triple) %in% germline$Suggested.Match)
res_germline <- subset(res_wt_triple, germlineEnriched)
ddsGermline <- dds1[rownames(dds1) %in% germline$Suggested.Match,]

# Extract upregulated and downregulated genes
up_wt_triple <- extract_genes(res_wt_triple, "up")
down_wt_triple <- extract_genes(res_wt_triple, "down")

# Extract upregulated and downregulated germline-enriched genes
up_wt_triple_germline <- extract_genes(res_germline, "up")
down_wt_triple_germline <- extract_genes(res_germline, "down")

# GO enrichment analysis for upregulated genes
genes_for_GO_up <- rownames(up_wt_triple)
perform_GO_enrichment(genes_for_GO_up)

# GO enrichment analysis for downregulated genes
genes_for_GO_down <- rownames(down_wt_triple)
perform_GO_enrichment(genes_for_GO_down)

# GO enrichment analysis for germline-enriched upregulated genes
genes_for_GO_up_germ <- rownames(up_wt_triple_germline)
perform_GO_enrichment(genes_for_GO_up_germ)

# GO enrichment analysis for germline-enriched downregulated genes
genes_for_GO_down_germ <- rownames(down_wt_triple_germline)
perform_GO_enrichment(genes_for_GO_down_germ)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Analysis of FBF2 bound ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create MA plot with FBF2 annotations
do_MA(res_wt_triple, FBF2bound, "mistyrose1", "darkorange") 

# Create heatmap for all genes and FBF2-bound genes
do_heatmap(res_wt_triple, dds1, coldata1, FBF2bound, "Is FBF2 bound", germline)

# Analyze upregulated genes and FBF2-bound genes
count_FBF_genes(res_wt_triple, FBF2bound, "up")
test_overlap(up_wt_triple, FBF2bound, dds1, "WT", "Triple", "FBF2 bound", "upregulated")

# Analyze downregulated genes and FBF2-bound genes
count_FBF_genes(res_wt_triple, FBF2bound, "down")
test_overlap(down_wt_triple, FBF2bound, dds1, "WT", "Triple", "FBF2 bound", "downregulated")

# ~~~~~~~~~~Germline enriched~~~~~~~~~~~

# Filter FBF2-bound genes for germline-enriched genes
FBF2Germline <- FBF2bound[FBF2bound$Suggested.Match %in% germline$Suggested.Match, ]

# Create heatmap for germline-enriched genes and FBF2-bound genes
do_heatmap(res_germline, ddsGermline, coldata1, FBF2bound, "Is FBF2 bound", germline)

# Analyze upregulated germline-enriched genes and FBF2-bound genes
count_FBF_genes(res_germline, FBF2bound, "up")
test_overlap(up_wt_triple_germline, FBF2Germline, ddsGermline, "WT", "Triple", "FBF2 bound in germline", "upregulated")

# Analyze downregulated germline-enriched genes and FBF2-bound genes
count_FBF_genes(res_germline, FBF2bound, "down")
test_overlap(down_wt_triple_germline, FBF2Germline, ddsGermline, "WT", "Triple", "FBF2 bound in germline", "downregulated")

# ~~~~~~~~~~GO~~~~~~~~~~~

# Extract FBF2-bound and upregulated genes
FBF2_up <- extract_FBF_genes(res_wt_triple, FBF2bound, "up")
# GO enrichment analysis for FBF2-bound and upregulated genes
print("GO enrichment for FBF2-bound and upregulated genes:")
perform_GO_enrichment(FBF2_up)

# Extract FBF2-bound and downregulated genes
FBF2_down <- extract_FBF_genes(res_wt_triple, FBF2bound, "down")
# GO enrichment analysis for FBF2-bound and downregulated genes
print("GO enrichment for FBF2-bound and downregulated genes:")
perform_GO_enrichment(FBF2_down)

# ~~~~~~~~~~analysis of only FBF2 bound for discussion~~~~~~~~~~~

# FBF2 bound and not FBF1 bound
FBF2_only <- FBF2bound[!(FBF2bound$Suggested.Match %in% FBF1bound$Suggested.Match), ]
FBF2_only_up <- extract_FBF_genes(res_wt_triple, FBF2_only, "up")
FBF2_only_down <- extract_FBF_genes(res_wt_triple, FBF2_only, "down")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Analysis of FBF1 bound ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create MA plot with FBF1 annotations
do_MA(res_wt_triple, FBF1bound, "mistyrose1", "darkorange", subName = "FBF1") 

# Create heatmap for all genes and FBF1-bound genes
do_heatmap(res_wt_triple, dds1, coldata1, FBF1bound, "Is FBF1 bound", germline)

# Analyze upregulated genes and FBF1-bound genes
count_FBF_genes(res_wt_triple, FBF1bound, "up")
test_overlap(up_wt_triple, FBF1bound, dds1, "WT", "Triple", "FBF1 bound", "upregulated")

# Analyze downregulated genes and FBF1-bound genes
count_FBF_genes(res_wt_triple, FBF1bound, "down")
test_overlap(down_wt_triple, FBF1bound, dds1, "WT", "Triple", "FBF1 bound", "downregulated")

# ~~~~~~~~~~Germline enriched~~~~~~~~~~~

# Filter FBF1-bound genes for germline-enriched genes
FBF1Germline <- FBF1bound[FBF1bound$Suggested.Match %in% germline$Suggested.Match, ]

# Create heatmap for germline-enriched genes and FBF1-bound genes
do_heatmap(res_germline, ddsGermline, coldata1, FBF1bound, "Is FBF1 bound", germline)

# Analyze upregulated germline-enriched genes and FBF1-bound genes
count_FBF_genes(res_germline, FBF1bound, "up")
test_overlap(up_wt_triple_germline, FBF1Germline, ddsGermline, "WT", "Triple", "FBF1 bound in germline", "upregulated")

# Analyze downregulated germline-enriched genes and FBF1-bound genes
count_FBF_genes(res_germline, FBF1bound, "down")
test_overlap(down_wt_triple_germline, FBF1Germline, ddsGermline, "WT", "Triple", "FBF1 bound in germline", "downregulated")

# ~~~~~~~~~~GO~~~~~~~~~~~

# Extract FBF1-bound and upregulated genes
FBF1_up <- extract_FBF_genes(res_wt_triple, FBF1bound, "up")
# GO enrichment analysis for FBF1-bound and upregulated genes
print("GO enrichment for FBF1-bound and upregulated genes:")
perform_GO_enrichment(FBF1_up)

# Extract FBF1-bound and downregulated genes
FBF1_down <- extract_FBF_genes(res_wt_triple, FBF1bound, "down")
# GO enrichment analysis for FBF1-bound and downregulated genes
print("GO enrichment for FBF1-bound and downregulated genes:")
perform_GO_enrichment(FBF1_down)

# ~~~~~~~~~~analysis of only FBF1 bound for discussion~~~~~~~~~~~

FBF1_only <- FBF1bound[!(FBF1bound$Suggested.Match %in% FBF2bound$Suggested.Match), ]
FBF1_only_up <- extract_FBF_genes(res_wt_triple, FBF1_only, "up")
FBF1_only_down <- extract_FBF_genes(res_wt_triple, FBF1_only, "down")

