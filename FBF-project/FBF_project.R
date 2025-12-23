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
#font_import()
#loadfonts(device = "win")
data(GeneOverlap)

setwd("C:\\Users\\gilad\\OneDrive\\Documents")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~constants~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
readsLim <- 10
lockBinding("readsLim", globalenv())

alphaLim <- 0.05
lockBinding("alphaLim", globalenv())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prep_counts <- function(file) {
  counts <- read.delim(file, header=TRUE, row.names = 1, sep = ",")
  row_names_to_remove <- c("WBGene00219609", "WBGene00219306", "WBGene00189951", "WBGene00001401") #triple, fbf1 "WBGene00001401"
  #row_names_to_remove <- c("WBGene00001402")
  counts <- counts[!(row.names(counts) %in% row_names_to_remove),]
  
  # counts <- counts[!grepl('?', counts$raw.WT_1, fixed = TRUE),]
  # names <- row.names(counts)
  # counts <- apply(counts,2,as.numeric)
  # row.names(counts) <- names
  
  return(counts[which(rowSums(counts) >= readsLim),])
}

prep_condition <- function() {
  #return(factor(c("Triple", "Triple", "Triple", "WT", "WT", "WT")))
  return(factor(c("Fbf1", "Fbf1", "Fbf1", "Fbf1_xtriple", "Fbf1_xtriple", "Fbf1_xtriple"))) 
  #return(factor(c("FBF2", "FBF2", "FBF2", "FBF2xTriple", "FBF2xTriple", "FBF2xTriple")))
  #return(factor(c("Fbf2_overexpression", "Fbf2_overexpression", "Fbf2_overexpression","WT", "WT", "WT", "WT")))
}

prep_coldata <- function(counts, condition) {
  return(data.frame(row.names = colnames(counts), condition))
}

#load data  
prepare_data <- function(file, refer) {
  counts <- prep_counts(file)
  condition <- prep_condition()
  coldata <- prep_coldata(counts, condition)
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
  dds$condition <- relevel(dds$condition, ref = refer)
  dds <- DESeq(dds)
  return(dds)
}  


#a modification of plotPCA function to match all PCs wanted
do_pca_helper <- function(object, first_pc, sec_pc, intgroup="condition", ntop=500, returnData=FALSE) {
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,first_pc], PC2=pca$x[,sec_pc], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(first_pc, sec_pc)] 
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC",first_pc," (",round(percentVar[first_pc] * 100),"%)")) +
    ylab(paste0("PC",sec_pc," (",round(percentVar[sec_pc] * 100),"%)")) +
                          theme(panel.background = element_blank(), 
                          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                          legend.title = element_blank(), text=element_text(size=12,  family="sans")) + 
                          coord_fixed(xlim=c(-15, 15), ylim = c(-15, 15)) +
                          ggforce::geom_mark_ellipse(aes(color = group, fill=group))+
                          scale_color_manual(values = c("blue","black")) +
                          scale_fill_manual(values = c("blue","black"))
  
}

#loads the PCA plot
do_pca <- function(dds, pc1, pc2) {
  transormed <- rlog(dds, blind=FALSE)
  do_pca_helper(transormed, pc1, pc2, intgroup="condition") 
}

#look at the results of dds
results_data <- function(dds, toCompare, baseLevel) {
  res <- results(dds, alpha = alphaLim, contrast = c("condition", toCompare, baseLevel)) 
  head(results(dds, tidy=TRUE))
  res <- res[order(res$padj),]
  head(res,10)
  summary(res)
  #write.csv(as.data.frame(res), file=glue('{toCompare}_vs_{baseLevel}_DGE.csv'))
  return(res)
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
  NSandFBF2bound = (data$name %in% subsetfile$Suggested.Match) & (data$sig == "NS")
  sigandFBF2bound = (data$name %in% subsetfile$Suggested.Match) & (data$sig != "NS")
  genes1 = subset(data, NSandFBF2bound)
  genes2 = subset(data, sigandFBF2bound)
  mean <- lfc <- sig <- name <- padj <-  NULL
  p <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig), size = size, alpha = alpha) +
    geom_point(data=genes1, color=subsetcolor1) +
    geom_point(data=genes2, color=subsetcolor2)
  
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
# Down up red
# FBF 2 bound blue
# checking what is the limit
# change the proportions - zoom in
#make MA plots for both shrunken and regular
do_MA <-function(res, subfile, subcolor1, subcolor2) {
  myggmaplot(res, subfile, subcolor1, subcolor2, main = expression("MA plot"),
           fdr = 0.05, fc = 0, size = 1, palette = c("mediumseagreen", "mediumseagreen", "darkgray"),
           genenames = as.vector(res$name),
           ggtheme = ggplot2::theme_minimal(),
           top = 0)
}

#make volcano plots, need to be done twice on shrunken and on regular
do_volcano <- function(res) {
  #Volcano Plot
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  resFilter<-
    with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-3,3)))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, padj<0.05 & abs(log2FoldChange)>0.5), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
  with(subset(res, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
  
  tmp <- res
  tmp$symbol <- rownames(tmp)
  res$symbol <- mapIds(org.Ce.eg.db, keys = rownames(res), keytype = "ENSEMBL", 
                           column = "SYMBOL")
  cond <- which(isNA(res$symbol))
  res[cond,]$symbol <- tmp[cond,]$symbol
  
  EnhancedVolcano(res,
                  lab = res$symbol,
                  x = 'log2FoldChange',
                  y = 'pvalue')
  
  EnhancedVolcano(res,
                  lab = res$symbol,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  boxedLabels = TRUE,
                  drawConnectors = TRUE,
                  widthConnectors = 0.3,
                  colConnectors = 'black',
                  max.overlaps = 50,
                  pCutoff = 10e-10,
                  FCcutoff = 0.5,
                  pointSize = 2.5,
                  labSize = 3.0)
  
}

#creating heat map
# orig arguments res, dds, coldata
do_heatmap <- function(res, dds, coldata, fileFBF, fileGermline) {

  sigs <- na.omit(res)
  sigs <- sigs[sigs$padj < alphaLim, ]
  sigs.df <- as.data.frame(sigs)
  # (sigs.df$baseMean > 50) 
  sigs.df <- sigs.df[(abs(sigs.df$log2FoldChange)> 0),]
  sigs.df <- sigs.df[order(sigs.df$log2FoldChange, decreasing = TRUE),]
  
  tmp <- sigs.df
  tmp$symbol <- rownames(tmp)
  sigs.df$symbol <- mapIds(org.Ce.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL", 
                           column = "SYMBOL")
  cond <- which(isNA(sigs.df$symbol))
  sigs.df[cond,]$symbol <- tmp[cond,]$symbol
  
  sigs.df$isGermlIne[rownames(sigs.df) %in% fileGermline$Suggested.Match] <- 1
  sigs.df$isGermlIne[!(rownames(sigs.df) %in% fileGermline$Suggested.Match)] <- 0
  sigs.df$isFBF2Bound[rownames(sigs.df) %in% fileFBF$Suggested.Match] <- 1
  sigs.df$isFBF2Bound[!(rownames(sigs.df) %in% fileFBF$Suggested.Match)] <- 0
  
  rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
  mat<-assay(rlog_out)[rownames(sigs.df), rownames(coldata)] #sig genes x samples
  colnames(mat) <- rownames(coldata)
  mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
  colnames(mat.scaled)<-colnames(mat)
  
  #tripleTmp <- subset(mat.scaled, select = c(Triple_1_S2, Triple_2_S9, Triple_3_S10))
  #tripleTmp <- as.data.frame(rowMeans(tripleTmp))
  #colnames(tripleTmp) <- "mean"
  #Triple <- tripleTmp$mean
  #mat.scaled <- cbind(mat.scaled, Triple)
  
  #wtTmp <- subset(mat.scaled, select = c(Wt_1_S1, WT_2_S7, WT_3_S8))
  #wtTmp <- as.data.frame(rowMeans(wtTmp))
  #colnames(wtTmp) <- "mean"
  #WT <- wtTmp$mean
  #mat.scaled <- cbind(mat.scaled, WT)

  #mat.scaled <- subset(mat.scaled, select = c(Triple, WT))
  
  
  isGermlIne <- as.matrix(sigs.df$isGermlIne)
  colnames(isGermlIne)<-"isGermlIne"
  isFBF2Bound <- as.matrix(sigs.df$isFBF2Bound)
  colnames(isFBF2Bound)<-"isFBF2Bound"
  
  colHeat <- colorRamp2(c(-1, 0, 1), c("blue4","white", "red3"))
  
  h1 <- Heatmap(mat.scaled, cluster_rows = T, cluster_columns = T, 
                column_labels = colnames(mat.scaled), name = "z_score", col = colHeat)
  
  h2 <- Heatmap(isGermlIne,
               cluster_rows = F, name="is Germline enriched", col=c("khaki", "forestgreen"))
  
  h3 <- Heatmap(isFBF2Bound,  
                cluster_rows = F, name = "is FBF2 Bound", col = c("gray", "maroon"))
  
  h<-h1+h2+h3
  #tiff("Heatmap_wt_FBF2_germline_cluster_by_samples_and genes.tiff", width = 1000, height = 1200, res = 150)
  print(h)
  #dev.off()
}

#count FBF and up
count_upregulated_FBF_genes <- function(res, FBfile) {
  res <- na.omit(res)
  up <- res[(res$log2FoldChange>0)&(res$padj<0.05),]
  
  FBFup <- sum(rownames(up) %in% FBfile$Suggested.Match)
  return(FBFup)
}

#count FBF and up
count_downregulated_FBF_genes <- function(res, FBfile) {
  res <- na.omit(res)
  down <- res[(res$log2FoldChange<0)&(res$padj<0.05),]
  
  FBFdown <- sum(rownames(down) %in% FBfile$Suggested.Match)
  return(FBFdown)
}

#extracting upregulated genes
extract_upregulated_genes <- function(res, refer, compareTo) {
  res <- na.omit(res)
  up <- res[(res$log2FoldChange>0)&(res$padj<0.05),]
  up$WB <- rownames(up)
  #write.csv(as.data.frame(up), file=glue('Upregulated_{refer}_vs_{compareTo}.csv'))
  return(up)
}

#extracting downregulated genes
extract_downregulated_genes <- function(res, refer, compareTo) {
  res <- na.omit(res)
  down <- res[(res$log2FoldChange<0)&(res$padj<0.05),]
  down$WB <- rownames(down)
  #write.csv(as.data.frame(down), file=glue('Downregulated_{refer}_vs_{compareTo}.csv'))
  return(down)
}

#extracting downregulated genes
extract_downregulated_all <- function(res) {
  res <- na.omit(res)
  down <- res[(res$log2FoldChange<0),]
  down$WB <- rownames(down)
  return(down)
}

#extracting downregulated genes
extract_downregulated_filtered <- function(res) {
  res <- na.omit(res)
  down <- res[(res$log2FoldChange<0) & (res$baseMean > 5) & (abs(res$log2FoldChange) > ((5/(res$baseMean^0.5))+0.3)),]
  down$WB <- rownames(down)
  return(down)
}

make_table <- function(res, dds, rawCounts, listForNames, FBFlist, germline) {
  table_deseq <- data.frame(row.names = listForNames$Gene_stable_ID)
  table_deseq$Gene_name <- listForNames$Gene_name[match(row.names(table_deseq), listForNames$Gene_stable_ID)]
  table_deseq$raw_FBF1_1 <- rawCounts$Fbf_1_S11[match(row.names(table_deseq), rawCounts$Geneid)]
  table_deseq$raw_FBF1_2 <- rawCounts$Fbf_1_S12[match(row.names(table_deseq), rawCounts$Geneid)]
  table_deseq$raw_FBF1_3 <- rawCounts$Fbf_1_S3[match(row.names(table_deseq), rawCounts$Geneid)]
  table_deseq$raw_FBF1Xtriple_1 <- rawCounts$Fbf_1xtriple_1_S2[match(row.names(table_deseq), rawCounts$Geneid)]
  table_deseq$raw_FBF1Xtriple_2 <- rawCounts$Fbf_1xtriple_2_S3[match(row.names(table_deseq), rawCounts$Geneid)]
  table_deseq$raw_FBF1Xtriple_3 <- rawCounts$Fbf_1xtriple_3_S4[match(row.names(table_deseq), rawCounts$Geneid)]
  
  normal <- estimateSizeFactors(dds)
  normalized_counts <- data.frame(counts(normal, normalized=TRUE))
  table_deseq$norm_FBF1_1 <- normalized_counts$Fbf_1_S11[match(row.names(table_deseq), row.names(normalized_counts))]
  table_deseq$norm_FBF1_2 <- normalized_counts$Fbf_1_S12[match(row.names(table_deseq), row.names(normalized_counts))]
  table_deseq$norm_FBF1_3 <- normalized_counts$Fbf_1_S3[match(row.names(table_deseq), row.names(normalized_counts))]
  table_deseq$norm_FBF1Xtriple_1 <- normalized_counts$Fbf_1xtriple_1_S2[match(row.names(table_deseq), row.names(normalized_counts))]
  table_deseq$norm_FBF1Xtriple_2 <- normalized_counts$Fbf_1xtriple_2_S3[match(row.names(table_deseq), row.names(normalized_counts))]
  table_deseq$norm_FBF1Xtriple_3 <- normalized_counts$Fbf_1xtriple_3_S4[match(row.names(table_deseq), row.names(normalized_counts))]
  
  table_deseq$baseMean <- res$baseMean[match(row.names(table_deseq), row.names(res))]
  table_deseq$log2FoldChange <- res$log2FoldChange[match(row.names(table_deseq), row.names(res))]
  table_deseq$lfcSE <- res$lfcSE[match(row.names(table_deseq), row.names(res))]
  table_deseq$stat <- res$stat[match(row.names(table_deseq), row.names(res))]
  table_deseq$pvalue <- res$pvalue[match(row.names(table_deseq), row.names(res))]
  table_deseq$padj <- res$padj[match(row.names(table_deseq), row.names(res))]
  
  table_deseq$Significantly_up[(table_deseq$padj < alphaLim) & (table_deseq$log2FoldChange > 0)] <- "True"
  table_deseq$Significantly_up[!((table_deseq$padj < alphaLim) & (table_deseq$log2FoldChange > 0))] <- "False"
  table_deseq$Significantly_down[(table_deseq$padj < alphaLim) & (table_deseq$log2FoldChange < 0)] <- "True"
  table_deseq$Significantly_down[!((table_deseq$padj < alphaLim) & (table_deseq$log2FoldChange < 0))] <- "False" 
  
  table_deseq$is_germlIne_enriched[rownames(table_deseq) %in% germline$Suggested.Match] <- "True"
  table_deseq$is_germlIne_enriched[!(rownames(table_deseq) %in% germline$Suggested.Match)] <- "False"
  table_deseq$is_FBF2_bound[rownames(table_deseq) %in% FBFlist$Suggested.Match] <- "True"
  table_deseq$is_FBF2_bound[!(rownames(table_deseq) %in% FBFlist$Suggested.Match)] <- "False"
  
  table_deseq[is.na(table_deseq)] <- "?"
  write.csv(as.data.frame(table_deseq), file=glue('FBF1_FBF1Xtriple_deseqRes.csv'))
  return(table_deseq)
}

#~~~~~~~~~~~~~prepare gene lists for hypergeometric test~~~~~~~~~~~
prepare_gene_lists <-function(file) {
  list <- read.delim(file, header=TRUE, row.names = NULL, sep = ",")
  colnames(list) <- colnames(list)[1:ncol(list)]
  list <- list[!grepl("not found", list$Suggested.Match),]
  list <- list[!(is.na(list$Suggested.Match) | list$Suggested.Match == ""),]
  return(list)
}

#preparing FBF1 bound list with WB names
prepare_FBF1_list <-function(file, file2) {
  list1 <- read.delim(file, header=TRUE, row.names = NULL, sep = ",")
  list2 <- read.delim(file2, header=TRUE, row.names = NULL, sep = ",")
  list1$Suggested.Match <- list2$Suggested.Match[match(list1$Gene.name, list2$Gene.name)]
  colnames(list1) <- colnames(list1)[1:ncol(list1)]
  list1 <- list1[!grepl("not found", list1$Suggested.Match),]
  list1 <- list1[!(is.na(list1$Suggested.Match) | list1$Suggested.Match == ""),]
  return(list1)
}

# create list of only fbf1 genes
create_only_fbf1 <- function(list1, list2) {
  only_FBF_1<- subset(list1,!(list1$Suggested.Match %in% list2$Suggested.Match))
  return(only_FBF_1)
}

# create lists of only fbf2 genes
create_only_fbf2 <- function(list1, list2) {
  only_FBF_2<-subset(list2,!(list2$Suggested.Match %in% list1$Suggested.Match))
  return(only_FBF_2)
}

# create lists of the intersection of fbf1 and fbf2
create_intersect_fbf1_fbf2 <- function(list1, list2) {
  intersect_FBF_1_and_FBF_2<-as.data.frame(intersect(list1$Suggested.Match,list2$Suggested.Match))
  colnames(intersect_FBF_1_and_FBF_2)<-"Suggested.Match"
  return(intersect_FBF_1_and_FBF_2)
}

# create lists of the intersection of fbf1 and fbf2
create_union_fbf1_fbf2 <- function(list1, list2) {
  both_FBF_1_and_FBF_2<-as.data.frame(union(list1$Suggested.Match,list2$Suggested.Match))
  colnames(both_FBF_1_and_FBF_2)<-"Suggested.Match"
  return(both_FBF_1_and_FBF_2)
}

#~~~~~~~~~~~~~~~~~hypergeometric test~~~~~~~~~~~~~~
test_overlap <- function(geneVecRegulated, FbfList, dds, refer, comp, fbf, regulation) {
  #construct a GeneOverlap object
  overlapObj <- newGeneOverlap(geneVecRegulated$WB,
                             FbfList$Suggested.Match,
                             genome.size = length(dds))
  
  #test overlap
  overlapObj <- testGeneOverlap(overlapObj)
  overlapObj
  print(overlapObj)
  overlapObj <- as.data.frame(getIntersection(overlapObj))
  
  #csv of overlap genes
  #write.csv(overlapObj, file=glue('intersection_of_{refer}_vs_{comp}_{regulation}_with_{fbf}.csv'))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~main~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FBF1bound <- prepare_FBF1_list("FBF-1_pull_down_RNA.csv", "FBF-1_conversion.csv")
FBF2bound <- prepare_gene_lists("FBF-2_pull_down_RNA.csv")
germline <-  read.delim("list_of_germline_enriched_genes_WB_ID _(reinke_list).csv", header=TRUE, row.names = NULL, sep = ",")
colnames(germline) <- c("Suggested.Match")
nameList <-  read.delim("yuval_analysis.csv", header=TRUE, row.names = NULL, sep = ",")

rawCounts <- read.delim("WT_vs_Triple_analysis.csv", header=TRUE, row.names = NULL, sep = ",")

only_FBF1 <- create_only_fbf1(FBF1bound, FBF2bound)
only_FBF2 <- create_only_fbf2(FBF1bound, FBF2bound)
intersect_FBF1_FBF2 <- create_intersect_fbf1_fbf2(FBF1bound, FBF2bound)
FBF1_FBF2_union <- create_union_fbf1_fbf2(FBF1bound, FBF2bound)


#~~~~~~~~~~~~~~WT_reference_for_05/03/23_analysis~~~~~~~~~~~~~~~~~~~
coldata1 <- prep_coldata(prep_counts("WT_vs_Triple_analysis.csv"), prep_condition())
dds1 <- prepare_data("WT_vs_Triple_analysis.csv", "WT")

do_pca(dds1, 1, 2)
plotDispEsts(dds1)

#~~~~~~~~~~~~~~~~~WT_VS_Triple~~~~~~~~~~~~~~~~~~~
res_wt_triple <- results_data(dds1, "Triple", "WT")
do_MA(res_wt_triple, FBF2bound, "mistyrose1", "darkorange") 
do_volcano(res_wt_triple)

germlineEnriched <- (rownames(res_wt_triple) %in% germline$Suggested.Match)
res_germline <- subset(res_wt_triple, germlineEnriched)
ddsGermline <- dds1[rownames(dds1) %in% germline$Suggested.Match,]

FBF2Germline <- FBF2bound$Suggested.Match %in% germline$Suggested.Match
FBF2Germline <- subset(FBF2bound, FBF2Germline)

tb <- make_table(res_wt_triple, dds1, rawCounts, nameList, FBF2bound, germline)

#all genes
do_heatmap(res_wt_triple, dds1, coldata1, FBF2bound, germline)

count_upregulated_FBF_genes(res_wt_triple, FBF2bound)
up_wt_triple <- extract_upregulated_genes(res_wt_triple, "WT", "Triple")
#up and Fbf2 bound
test_overlap(up_wt_triple, FBF2bound, dds1, "WT", "Triple", "fbf2_bound", "upregulated")

count_downregulated_FBF_genes(res_wt_triple, FBF2bound)
down_wt_triple <- extract_downregulated_genes(res_wt_triple, "WT", "Triple")
#down and Fbf2 bound
test_overlap(down_wt_triple, FBF2bound, dds1, "WT", "Triple", "fbf2_bound", "downregulated")

#germline enriched
do_heatmap(res_germline, ddsGermline, coldata1, FBF2bound, germline)
count_upregulated_FBF_genes(res_germline, FBF2bound)
up_wt_triple_germline <- extract_upregulated_genes(res_germline, "WT", "Triple")
#up and Fbf2 bound of germline enriched
test_overlap(up_wt_triple_germline, FBF2bound, ddsGermline, "WT", "Triple", "fbf2_bound", "upregulated")

count_downregulated_FBF_genes(res_germline, FBF2bound)
down_wt_triple_germline <- extract_downregulated_genes(res_germline, "WT", "Triple")
#down and Fbf2 bound of germline enriched
test_overlap(down_wt_triple_germline, FBF2bound, ddsGermline, "WT", "Triple", "fbf2_bound", "downregulated")


#germline enriched and FBF2 in germline
#up and Fbf2 bound of germline enriched and FBF2 in germline
test_overlap(up_wt_triple_germline, FBF2Germline, ddsGermline, "WT", "Triple", "fbf2_bound", "upregulated")

#down and Fbf2 bound of germline enriched and FBF2 in germline
test_overlap(down_wt_triple_germline, FBF2Germline, ddsGermline, "WT", "Triple", "fbf2_bound", "downregulated")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#up and Fbf1 bound
test_overlap(up_wt_triple, FBF1bound, dds1, "WT", "Triple", "fbf1_bound", "upregulated")
#up and Fbf1 specific
test_overlap(up_wt_triple, only_FBF1, dds1, "WT", "Triple", "fbf1_specific", "upregulated")

#up and Fbf2 specific
test_overlap(up_wt_triple, only_FBF2, dds1, "WT", "Triple", "fbf2_specific", "upregulated")
#up and fbf1 fbf2 intersection 
test_overlap(up_wt_triple, intersect_FBF1_FBF2, dds1, "WT", "Triple", "fbf1_fbf2_intersection", "upregulated")
#up and fbf1 fbf2 union 
test_overlap(up_wt_triple, FBF1_FBF2_union, dds1, "WT", "Triple", "fbf1_fbf2_union", "upregulated")

#down and Fbf1 bound
test_overlap(down_wt_triple, FBF1bound, dds1, "WT", "Triple", "fbf1_bound", "downregulated")
#down and Fbf1 specific
test_overlap(down_wt_triple, only_FBF1, dds1, "WT", "Triple", "fbf1_specific", "downregulated")
#down and Fbf2 specific
test_overlap(down_wt_triple, only_FBF2, dds1, "WT", "Triple", "fbf2_specific", "downregulated")
#down and fbf1 fbf2 intersection 
test_overlap(down_wt_triple, intersect_FBF1_FBF2, dds1, "WT", "Triple", "fbf1_fbf2_intersection", "downregulated")
#down and fbf1 fbf2 union 
test_overlap(down_wt_triple, FBF1_FBF2_union, dds1, "WT", "Triple", "fbf1_fbf2_union", "downregulated")


#~~~~~~~~~~~~~~WT_reference_for_23/03/23_analysis~~~~~~~~~~~~~~~~~~~
coldata2 <- prep_coldata(prep_counts("WT_vs_FBF2_analysis.csv"), prep_condition())
dds2 <- prepare_data("WT_vs_FBF2_analysis.csv", "WT")

do_pca(dds2, 1, 2)

#~~~~~~~~~~~~~~~~~WT_VS_FBF2~~~~~~~~~~~~~~~~~~~
res_wt_FBF2 <- results_data(dds2, "Fbf2", "WT")

germlineEnrichedFBF2 <- (rownames(res_wt_FBF2) %in% germline$Suggested.Match)
res_FBF2_germline <- subset(res_wt_FBF2, germlineEnrichedFBF2)
dds2Germline <- dds2[rownames(dds2) %in% germline$Suggested.Match,]

rawCountsFBF2 <- read.delim("WT_vs_FBF2_analysis.csv", header=TRUE, row.names = NULL, sep = ",")
tbFBF2 <- make_table(res_wt_FBF2, dds2, rawCountsFBF2, nameList, FBF2bound, germline)

#all genes
do_heatmap(res_wt_FBF2, dds2, coldata2, FBF2bound, germline)
#germline
do_heatmap(res_FBF2_germline, dds2Germline, coldata2, FBF2bound, germline)

count_upregulated_FBF_genes(res_wt_FBF2, FBF2bound)
up_wt_FBF2 <- extract_upregulated_genes(res_wt_FBF2, "WT", "FBF2")
#up and Fbf2 bound
test_overlap(up_wt_FBF2, FBF2bound, dds2, "WT", "FBF2", "fbf2_bound", "upregulated")

count_downregulated_FBF_genes(res_wt_FBF2, FBF2bound)
down_wt_FBF2 <- extract_downregulated_genes(res_wt_FBF2, "WT", "FBF2")
#down and Fbf2 bound
test_overlap(down_wt_FBF2, FBF2bound, dds2, "WT", "FBF2", "fbf2_bound", "downregulated")



#~~~~~~~~~~~~~~FBF2_reference_for_30/07/23_analysis~~~~~~~~~~~~~~~~~~~
coldata3 <- prep_coldata(prep_counts("FBF2_VS_FBF2xTriple_analysis.csv"), prep_condition())
dds3 <- prepare_data("FBF2_VS_FBF2xTriple_analysis.csv", "FBF2")
do_pca(dds3, 1, 2)

#~~~~~~~~~~~~~~~~~FBF2_VS_FBF2xTriple~~~~~~~~~~~~~~~~~~~
res_FBF2 <- results_data(dds3, "FBF2xTriple", "FBF2")

germlineEnriched3 <- (rownames(res_FBF2) %in% germline$Suggested.Match)
res_germline3 <- subset(res_FBF2, germlineEnriched3)
ddsGermline3 <- dds3[rownames(dds3) %in% germline$Suggested.Match,]

FBF1Germline <- FBF1bound$Suggested.Match %in% germline$Suggested.Match
FBF1Germline <- subset(FBF1bound, FBF1Germline)

#all genes
#do_heatmap(res_FBF2, dds3, coldata3, FBF1bound, germline)

count_upregulated_FBF_genes(res_FBF2, FBF1bound)
up_FBF2_FBF2xTriple <- extract_upregulated_genes(res_FBF2, "FBF2", "FBF2xTriple")
#up and Fbf1 bound
test_overlap(up_FBF2_FBF2xTriple, FBF1bound, dds3, "FBF2", "FBF2xTriple", "fbf1_bound", "upregulated")

count_downregulated_FBF_genes(res_FBF2, FBF1bound)
down_FBF2_FBF2xTriple <- extract_downregulated_genes(res_FBF2, "FBF2", "FBF2xTriple")
#down and Fbf1 bound
test_overlap(down_FBF2_FBF2xTriple, FBF1bound, dds3, "FBF2", "FBF2xTriple", "fbf1_bound", "downregulated")


#germline enriched
#do_heatmap(res_germline3, ddsGermline3, coldata3, FBF1bound, germline)
count_upregulated_FBF_genes(res_germline3, FBF1bound)
up_FBF2_FBF2xTriple_germline <- extract_upregulated_genes(res_germline3, "FBF2", "FBF2xTriple")
#up and Fbf1 bound of germline enriched and FBF1 in germline
test_overlap(up_FBF2_FBF2xTriple_germline, FBF1Germline, ddsGermline3, "FBF2", "FBF2xTriple", "fbf1_bound", "upregulated")

count_downregulated_FBF_genes(res_germline3, FBF1bound)
down_FBF2_FBF2xTriple_germline <- extract_downregulated_genes(res_germline3, "FBF2", "FBF2xTriple")
#down and Fbf1 bound of germline enriched and FBF1 in germline
test_overlap(down_FBF2_FBF2xTriple_germline, FBF1Germline, ddsGermline3, "FBF2", "FBF2xTriple", "fbf1_bound", "downregulated")


#~~~~~~~~~~~~~~FBF1_reference_for_12/05/24_analysis~~~~~~~~~~~~~~~~~~~
coldata3 <- prep_coldata(prep_counts("fbf1_vs_fbf1triple.csv"), prep_condition())
dds3 <- prepare_data("fbf1_vs_fbf1triple.csv", "Fbf1")
do_pca(dds3, 1, 2)


#~~~~~~~~~~~~~~~~~FBF1_VS_FBF1xTriple~~~~~~~~~~~~~~~~~~~
res_FBF1_FBf1Triple <- results_data(dds3, "Fbf1_xtriple", "Fbf1")

germlineEnriched3 <- (rownames(res_FBF1_FBf1Triple) %in% germline$Suggested.Match)
res_germline3 <- subset(res_FBF1_FBf1Triple, germlineEnriched3)
ddsGermline3 <- dds3[rownames(dds3) %in% germline$Suggested.Match,]

FBF2Germline <- FBF2bound$Suggested.Match %in% germline$Suggested.Match
FBF2Germline <- subset(FBF2bound, FBF2Germline)

rawCounts_FBF1_FBF1Xtriple <- read.delim("fbf1_vs_fbf1triple.csv", header=TRUE, row.names = NULL, sep = ",")
tb3 <- make_table(res_FBF1_FBf1Triple, dds3, rawCounts_FBF1_FBF1Xtriple, nameList, FBF2bound, germline)

#all genes
do_heatmap(res_FBF1_FBf1Triple, dds3, coldata3, FBF2bound, germline)
count_upregulated_FBF_genes(res_FBF1_FBf1Triple, FBF2bound)
up_FBF1_FBF1xTriple <- extract_upregulated_genes(res_FBF1_FBf1Triple, "Fbf1", "Fbf1_xtriple")
#up and Fbf2 bound
test_overlap(up_FBF1_FBF1xTriple, FBF2bound, dds3, "Fbf1", "Fbf1_xtriple", "fbf2_bound", "upregulated")

count_downregulated_FBF_genes(res_FBF1_FBf1Triple, FBF2bound)
down_FBF1_FBF1xTriple <- extract_downregulated_genes(res_FBF1_FBf1Triple, "Fbf1", "Fbf1_xtriple")
#down and Fbf2 bound
test_overlap(down_FBF1_FBF1xTriple, FBF2bound, dds3, "Fbf1", "Fbf1_xtriple", "fbf2_bound", "downregulated")


#germline enriched
#do_heatmap(res_germline3, ddsGermline3, coldata3, FBF1bound, germline)
count_upregulated_FBF_genes(res_germline3, FBF2bound)
up_FBF1_FBF1xTriple_germline <- extract_upregulated_genes(res_germline3, "Fbf1", "Fbf1_xtriple")
#up and Fbf1 bound of germline enriched and FBF1 in germline
test_overlap(up_FBF1_FBF1xTriple_germline, FBF2Germline, ddsGermline3, "Fbf1", "Fbf1_xtriple", "fbf2_bound", "upregulated")

count_downregulated_FBF_genes(res_germline3, FBF2bound)
down_FBF1_FBF1xTriple_germline <- extract_downregulated_genes(res_germline3,"Fbf1", "Fbf1_xtriple")
#down and Fbf1 bound of germline enriched and FBF1 in germline
test_overlap(down_FBF1_FBF1xTriple_germline, FBF2Germline, ddsGermline3, "Fbf1", "Fbf1_xtriple", "fbf2_bound", "downregulated")









#~~~~~~~~~~~~~~~~~~~~~~~FBF2 OVEREXPRESSION VS WT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rawCountsFBF2OE <- read.delim("FBF2_oe_VS_WT.csv", header=TRUE, row.names = NULL, sep = ",")

#~~~~~~~~~~~~~~WT_reference_for_17/10/23_analysis~~~~~~~~~~~~~~~~~~~
coldata1 <- prep_coldata(prep_counts("FBF2_oe_VS_WT.csv"), prep_condition())
dds1 <- prepare_data("FBF2_oe_VS_WT.csv", "WT")

do_pca(dds1, 1, 2)

res_wt_FBF2_OE <- results_data(dds1, "Fbf2_overexpression", "WT")

#tb <- make_table(res_wt_FBF2_OE, dds1, rawCountsFBF2OE, nameList, FBF2bound, germline)
down_wt_FBF2OE <- extract_downregulated_genes(res_wt_FBF2_OE, "WT", "FBF2OE")

down_wt_FBF2OE_all <- extract_downregulated_all(res_wt_FBF2_OE)

down_wt_FBF2OE_filtered <- extract_downregulated_filtered(res_wt_FBF2_OE)

#~~~~~~~~~~ WT Triple~~~~~~~~~~~
coldata2 <- prep_coldata(prep_counts("WT_vs_Triple_analysis.csv"), prep_condition())
dds2 <- prepare_data("WT_vs_Triple_analysis.csv", "WT")

res_wt_triple <- results_data(dds2, "Triple", "WT")

down_wt_triple <- extract_downregulated_genes(res_wt_triple, "WT", "Triple")

down_wt_triple_all <- extract_downregulated_all(res_wt_triple)

down_wt_triple_filtered <- extract_downregulated_filtered(res_wt_triple)

#~~~~~~~~~~~~~~~down triple and down fbf2 overexpression~~~~~~~~~~~~
overlapObj <- newGeneOverlap(down_wt_triple_filtered$WB,
                             down_wt_FBF2OE_filtered$WB,
                             genome.size = length(dds1))

#test overlap
overlapObj <- testGeneOverlap(overlapObj)
overlapObj
print(overlapObj)
overlapObj <- as.data.frame(getIntersection(overlapObj))
overlapObj$Gene_name <- nameList$Gene_name[match(overlapObj$'getIntersection(overlapObj)', nameList$Gene_stable_ID)]

#csv of overlap genes
write.csv(overlapObj, file='FBF2 OE downregulated overlapped with Triple downregulated filtered.csv')


