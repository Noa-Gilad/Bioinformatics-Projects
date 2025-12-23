library(DESeq2)
library(ggplot2)
library(glue)
library(EnhancedVolcano)
library(org.Hs.eg.db)
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
library(AnnotationDbi)
library(clusterProfiler)
library(msigdbr)

#font_import()
#loadfonts(device = "win")
data(GeneOverlap)

library("sva")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~constants~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
readsLim <- 300
lockBinding("readsLim", globalenv())

alphaLim <- 0.05
lockBinding("alphaLim", globalenv())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prep_counts <- function(file, toRemove) {
  counts <- read.delim(file, header=TRUE, row.names = 1, sep = ",")
  counts <- counts[!(row.names(counts) %in% row.names(toRemove)),] #to change
  return(counts[which(rowSums(counts) >= readsLim),])
}


prep_coldata <- function(file) {
  design_file = read.csv(file, sep=',', row.names = 1, stringsAsFactors = F)
  return(as.data.frame(apply(design_file, 2, factor)))
}


#load data
prepare_data <- function(file, refer, design_f, toRemove) {
  counts <- prep_counts(file, toRemove)
  coldata <- prep_coldata(design_f)
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition + sex)
  dds$condition <- relevel(dds$condition, ref = refer)
  dds <- DESeq(dds)
  return(dds)
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
  data$name_prefix <- sub("\\..*", "", data$name) # for this project only
  
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
  NSandInSubfile = (data$name_prefix %in% subsetfile$Suggested.Match) & (data$sig == "NS") #to change
  sigandinsubfile = (data$name_prefix %in% subsetfile$Suggested.Match) & (data$sig != "NS")
  genes1 = subset(data, NSandInSubfile)
  genes2 = subset(data, sigandinsubfile)
  
  tmp <- genes2
  tmp$symbol <- tmp$name_prefix
  genes2$symbol <- mapIds(org.Hs.eg.db, keys = genes2$name_prefix, keytype = "ENSEMBL", 
                       column = "SYMBOL")
  cond <- which(isNA(genes2$symbol))
  genes2[cond,]$symbol <- tmp[cond,]$symbol
  
  mean <- lfc <- sig <- name <- padj <-  NULL
  p <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig), size = size, alpha = alpha) +
    geom_point(data=genes1, color=subsetcolor1) + 
    geom_point(data=genes2, color=subsetcolor2, size=2.5) + 
    geom_label_repel(data = genes2, aes(label = symbol), 
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

#a modification of plotPCA function to match all PCs wanted
do_pca_helper <- function(object, first_pc, sec_pc, intgroup="condition", ntop=500, returnData=FALSE, 
                          pca_colors = c("black","blue")) {
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
    ggforce::geom_mark_ellipse(aes(color = group, fill=group)) +
    scale_color_manual(values = pca_colors) +
    scale_fill_manual(values = pca_colors) +
    geom_text(data = d, aes(label = name), vjust = -0.5, hjust = 1)
  
}

#loads the PCA plot
do_pca <- function(dds, pc1, pc2, p_colors, intg) {
  transormed <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5))
  do_pca_helper(transormed, pc1, pc2, intgroup=intg, pca_colors = p_colors) 
}

#creating heat map
# orig arguments res, dds, coldata
do_heatmap <- function(res, dds, coldata, subfile) {
  
  sigs <- na.omit(res)
  sigs <- sigs[sigs$padj < alphaLim, ]
  sigs.df <- as.data.frame(sigs)
  sigs.df <- sigs.df[(abs(sigs.df$log2FoldChange)> 1),]
  sigs.df <- sigs.df[order(sigs.df$log2FoldChange, decreasing = TRUE),]
  sigs.df$name_prefix <- sub("\\..*", "", rownames(sigs.df))
  
  sigs.df$isInsubfile[sigs.df$name_prefix %in% subfile$Suggested.Match] <- 1
  sigs.df$isInsubfile[!(sigs.df$name_prefix %in% subfile$Suggested.Match)] <- 0
  
  transformed <- vst(dds, blind=FALSE, nsub=sum(rowMeans(counts(dds, normalized=TRUE)) > 5)) #get normalized count data from dds object
  mat<-assay(transformed)[rownames(sigs.df), rownames(coldata)] #sig genes x samples
  colnames(mat) <- rownames(coldata)
  mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
  colnames(mat.scaled)<-colnames(mat)

  isInsubfile <- as.matrix(sigs.df$isInsubfile)
  colnames(isInsubfile)<-"isInsubfile"
  
  colHeat <- colorRamp2(c(-1, 0, 1), c("blue4","white", "red3"))
  
  h1 <- Heatmap(mat.scaled, cluster_rows = F, cluster_columns = T, 
                column_labels = colnames(mat.scaled), name = "z_score", col = colHeat, show_row_names = FALSE)
  
  h3 <- Heatmap(isInsubfile,  
                cluster_rows = F, name = "known biomarker", col = c("gray", "maroon"))
  
  h<-h1+h3

  print(h)

}

#extracting upregulated genes
extract_upregulated_genes <- function(res, refer, compareTo) {
  res <- na.omit(res)
  up <- res[(res$log2FoldChange>1)&(res$padj<0.05),]
  
  up$name_prefix <- sub("\\..*", "", rownames(up))
  tmp <- up
  tmp$symbol <- tmp$name_prefix
  up$symbol <- mapIds(org.Hs.eg.db, keys = up$name_prefix, keytype = "ENSEMBL", 
                          column = "SYMBOL")
  cond <- which(isNA(up$symbol))
  up[cond,]$symbol <- tmp[cond,]$symbol
  
  # write.csv(as.data.frame(up), file=glue('Upregulated_{refer}_vs_{compareTo}.csv'))
  return(up)
}

#extracting upregulated genes
extract_downregulated_genes <- function(res, refer, compareTo) {
  res <- na.omit(res)
  down <- res[(res$log2FoldChange< -1)&(res$padj<0.05),]
  
  down$name_prefix <- sub("\\..*", "", rownames(down))
  tmp <- down
  tmp$symbol <- tmp$name_prefix
  down$symbol <- mapIds(org.Hs.eg.db, keys = down$name_prefix, keytype = "ENSEMBL", 
                      column = "SYMBOL")
  cond <- which(isNA(down$symbol))
  down[cond,]$symbol <- tmp[cond,]$symbol
  
  #write.csv(as.data.frame(down), file=glue('Downregulated_{refer}_vs_{compareTo}.csv'))
  return(down)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~main~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colData <- prep_coldata( "design.csv")

remove <- read.delim("hemoglobin_genes.csv", header=TRUE, row.names = 1, sep = ",")
dds <- prepare_data("mRNA_count_table_all.csv", "control", "design.csv", remove)
do_pca(dds, 1, 2, c("black","blue"), c("sex"))

dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition + sex, colData(dds))
mod0 <- model.matrix(~ sex, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 18)

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
ddssva$SV4 <- svseq$sv[,4]
ddssva$SV5 <- svseq$sv[,5]
ddssva$SV6 <- svseq$sv[,6]
ddssva$SV7 <- svseq$sv[,7]
ddssva$SV8 <- svseq$sv[,8]
ddssva$SV9 <- svseq$sv[,9]
ddssva$SV10 <- svseq$sv[,10]
ddssva$SV11 <- svseq$sv[,11]
ddssva$SV12 <- svseq$sv[,12]
ddssva$SV13 <- svseq$sv[,13]
ddssva$SV14 <- svseq$sv[,14]
ddssva$SV15 <- svseq$sv[,15]
ddssva$SV16 <- svseq$sv[,16]
ddssva$SV17 <- svseq$sv[,17]
ddssva$SV18 <- svseq$sv[,18]
design(ddssva) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + 
  SV11 + SV12 + SV13 + SV14 + SV15 + SV16 + SV17 + SV18 + condition + sex
ddssva <- DESeq(ddssva)

do_pca(ddssva, 1, 2, c("black","blue"), c("condition"))

res <- results_data(ddssva, "patient", "control")

biomarkers <- read.delim("biomarkers.csv", sep = ",")

do_MA(res, biomarkers, "#1dce00", "#ec7e0a")
do_heatmap(res, ddssva, colData, biomarkers) # foldchange > 1

up <- extract_upregulated_genes(res, "patient", "control")
down <- extract_downregulated_genes(res, "patient", "control")

ddsvaInBio <- ddssva[sub("\\..*", "", rownames(ddssva)) %in% biomarkers$Suggested.Match,]
do_pca(ddsvaInBio, 1, 2, c("black","blue"), c("condition"))


ddsvaInBioUpDown <- ddsvaInBio[(rownames(ddsvaInBio) %in% rownames(up)) | (rownames(ddsvaInBio) %in%  rownames(down)), ]
do_pca(ddsvaInBioUpDown, 1, 2, c("black","blue"), c("condition"))

ddsvaInUpDown <- ddssva[(rownames(ddssva) %in% rownames(up)) | (rownames(ddssva) %in%  rownames(down)), ]
do_pca(ddsvaInUpDown, 1, 2, c("black","blue"), c("condition"))

#go analysis
sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05 & sigs$baseMean > 50,]
rownames(sigs) <- sub("\\..*", "", rownames(sigs))
genes_to_test_up <- rownames(sigs[sigs$log2FoldChange > 0,])
genes_to_test_down <- rownames(sigs[sigs$log2FoldChange < 0,])

GO_results_BP_up <- enrichGO(gene = genes_to_test_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results_BP_up)
fit_up <- plot(barplot(GO_results_BP_up, showCategory = 15))
png("BP_up.png", res = 250, width = 1400, height = 1800)
 
GO_results_BP_down <- enrichGO(gene = genes_to_test_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results_BP_down)
fit_down <- plot(barplot(GO_results_BP_down, showCategory = 15))
png("BP_down.png", res = 250, width = 1400, height = 1800)

#msigdb

cancer_gene_sets = msigdbr(species = "human", category = "C4")
msigdbr_t2g = cancer_gene_sets %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()

msigdb_up <- enricher(gene = genes_to_test_up, TERM2GENE = msigdbr_t2g)
as.data.frame(msigdb_up)
fit_msigdb_up <- plot(barplot(msigdb_up, showCategory = 20))

msigdb_down <- enricher(gene = genes_to_test_down, TERM2GENE = msigdbr_t2g)
as.data.frame(msigdb_down)
fit_msigdb_down <- plot(barplot(msigdb_down, showCategory = 20))

 
#C6

cancer_gene_sets_C6 = msigdbr(species = "human", category = "C6")
msigdbr_t2g_onco = cancer_gene_sets_C6 %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()

msigdb_up_C6 <- enricher(gene = genes_to_test_up, TERM2GENE = msigdbr_t2g_onco)
as.data.frame(msigdb_up_C6)
fit_msigdb_C6_up <- plot(barplot(msigdb_up_C6, showCategory = 20))

msigdb_down_C6 <- enricher(gene = genes_to_test_down, TERM2GENE = msigdbr_t2g_onco)
as.data.frame(msigdb_down_C6)
fit_msigdb_C6_up <- plot(barplot(msigdb_down_C6, showCategory = 20))

#C2

cancer_gene_sets_C2 = msigdbr(species = "human", category = "C2")
msigdbr_t2g_c2 = cancer_gene_sets_C2 %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()

msigdb_up_C2 <- enricher(gene = genes_to_test_up, TERM2GENE = msigdbr_t2g_c2)
as.data.frame(msigdb_up_C2)
fit_msigdb_C2_up <- plot(barplot(msigdb_up_C2, showCategory = 18))

msigdb_down_C2 <- enricher(gene = genes_to_test_down, TERM2GENE = msigdbr_t2g_c2)
as.data.frame(msigdb_down_C2)
fit_msigdb_C2_up <- plot(barplot(msigdb_down_C2, showCategory = 15))

#H

cancer_gene_sets_H = msigdbr(species = "human", category = "H")
msigdbr_t2g_H = cancer_gene_sets_H %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()

msigdb_up_H <- enricher(gene = genes_to_test_up, TERM2GENE = msigdbr_t2g_H)
as.data.frame(msigdb_up_H)
fit_msigdb_H_up <- plot(barplot(msigdb_up_H, showCategory = 20))

msigdb_down_H <- enricher(gene = genes_to_test_down, TERM2GENE = msigdbr_t2g_H)
as.data.frame(msigdb_down_H)
fit_msigdb_H_up <- plot(barplot(msigdb_down_H, showCategory = 20))

