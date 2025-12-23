###############################################################################
# Oocyte Aging Analysis – 12 Oocytes Dataset (Yuan et al.)
#
# Author: Noa Gilad
#
# Description:
#   Differential expression and quality analysis of human oocytes
#   comparing young vs old donors, based on the 12-oocyte dataset
#   reported in Yuan et al.
#
#   The script:
#     - Reproduces analysis using the original paper count table
#     - Performs an independent re-analysis using in-house generated counts
#     - Conducts DESeq2 normalization and differential expression
#     - Identifies and removes outlier oocytes (as reported in the paper)
#     - Generates PCA, MA plots, and expression QC visualizations
#     - Extracts age-associated DE genes and performs GO enrichment
#
# Input:
#   - Paper count table (paper_counts_12.txt)
#   - In-house count table (count_12_FC.txt)
#   - Sample metadata CSVs (paper_12_des.csv, design_12.csv)
#
# Output:
#   - PCA and QC plots
#   - DESeq2 result tables (young vs old)
#   - Gene lists for downstream enrichment analyses
#
# Notes:
#   - This script sources shared DESeq2 and QC helper pipelines.
#   - Outliers are removed to match the original publication.
#   - Designed for reproducibility and comparison with published results.
#
###############################################################################

# ~~~~~~~~~~~~ source pipeline ~~~~~~~~~~
source( "C:\\Users\\gilad\\R\\pipelines\\deseq2_pipline.R", local = FALSE)
source("C:\\Users\\gilad\\R\\pipelines\\preprocessing_and_QA.R", local = FALSE)

#  ~~~~~~~~~~~~~~~~~~~~~~~~~ analyse 12 oocytes ~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~ Testing paper's count table ~~~~~~~~~~~~
df_12_paper <- read_tsv("paper_counts_12.txt")
df_12_clean_pap <- df_12_paper %>% filter(!grepl("^__", gene_id))
df_12_pap_clean_with_rownames <- df_12_clean_pap %>% column_to_rownames(var = "gene_id")
colData_12_pap <- prep_coldata("paper_12_des.csv", condition_col = "age", ref_level = "Y")

counts_12_pap <- prep_counts(counts = df_12_pap_clean_with_rownames, readsLim = 10)
counts_12_pap <- counts_12_pap[, rownames(colData_12_pap)]
dds_12_pap <- prepare_data(counts = counts_12_pap, coldata = colData_12_pap, refer = "Y", des = ~ age, condition = "age")

transformed_12_pap <-  vst(dds_12_pap, blind=FALSE)
do_pca(object = transformed_12_pap, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"))

outlier_samples_12 <- c("YOF", "YOG", "YOB")

# Remove the outlier from the count table
counts_12_pap <- counts_12_pap[, !(colnames(counts_12_pap) %in% outlier_samples_12)]

# Remove the outlier from the colData
colData_12_pap <- colData_12_pap[!(rownames(colData_12_pap) %in% outlier_samples_12), , drop = FALSE]

dds_12_pap <- prepare_data(counts = counts_12_pap, coldata = colData_12_pap, refer = "Y", des = ~ age, condition = "age")

transformed_12_pap <-  vst(dds_12_pap, blind=FALSE)

do_pca(object = transformed_12_pap, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"))

res_12_pap <- results_data(dds = dds_12_pap, toCompare = "O", baseLevel = "Y", lfc = 0, alphaLim = 0.05, condition = "age", writeFile = FALSE)

extract_DEG_deseq2_data(res = res_12_pap, direction = "up", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_up_12_paper_deseq2_data")
extract_DEG_deseq2_data(res = res_12_pap, direction = "down", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_down_12_paper_deseq2_data")

do_MA_plot(res = res_12_pap, gene_list = gene_list, fc = 2, fdr = 0.05)


# ~~~~~~~~~~~~ Count table created in-house ~~~~~~~~~~~~~~

# df_12_oocytes <- read_tsv("count_matrix_12_oocytes.txt")
df_12_oocytes <- read_tsv("count_12_FC.txt", comment = "#")
df_12_clean <- df_12_oocytes %>% dplyr::select(Geneid, starts_with("SRR"))
colnames(df_12_clean) <- sub("Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(df_12_clean))
#df_12_clean <- df_12_oocytes %>% filter(!grepl("^__", gene_id))
#df_12_clean_with_rownames <- df_12_clean %>% column_to_rownames(var = "gene_id")
df_12_clean_with_rownames <- df_12_clean %>% column_to_rownames(var = "Geneid")

colData_12_exp <- prep_coldata("design_12.csv", condition_col = "age", ref_level = "Y")

counts_12 <- prep_counts(counts = df_12_clean_with_rownames, readsLim = 20)
counts_12 <- counts_12[, rownames(colData_12_exp)]
create_count_histogram(count_table = counts_12)

dds_12 <- DESeqDataSetFromMatrix(countData = counts_12, colData = colData_12_exp, design = ~ age)
dds_12 <- estimateSizeFactors(dds_12)
normalized_counts_12 <- counts(dds_12, normalized = TRUE)
log_transformed_12 <- as.data.frame(log2(normalized_counts_12 + 1))

boxplot(log_transformed_12, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,           # Make sample labels vertical
        col="lightblue", # Box color
        border="darkblue",
        main="RNA-Seq Expression Data",
        cex.axis=0.7,    # Reduce axis label size if needed
        par(mar=c(10,4,4,2) + 0.1)) # Increase bottom margin for labels

plot_all_pairs(log_counts = log_transformed_12, coldata = colData_12_exp, out_dir = "C:\\Users\\gilad\\R\\oocytes_aging\\12_results\\pairwise")

dds_12 <- prepare_data(counts = counts_12, coldata = colData_12_exp, refer = "Y", des = ~ age, condition = "age")

transformed_12 <-  vst(dds_12, blind=FALSE)
plot_vst_boxplot(transformed_12)

do_pca(dds = dds_12, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"))

# removing SRR12330964, SRR12330963, SRR12330965 like in the original paper 
outlier_samples_12 <- c("SRR12330964", "SRR12330963", "SRR12330965")

# Remove the outlier from the count table
df_12_clean_with_rownames <- df_12_clean_with_rownames[, !(colnames(df_12_clean_with_rownames) %in% outlier_samples_12)]
counts_12 <- prep_counts(counts = df_12_clean_with_rownames, readsLim = 10)

# Remove the outlier from colData
colData_12_exp <- colData_12_exp[!(rownames(colData_12_exp) %in% outlier_samples_12), , drop = FALSE]

dds_12 <- prepare_data(counts = counts_12, coldata = colData_12_exp, refer = "Y", des = ~ age, condition = "age")

transformed_12 <-  vst(dds_12, blind=FALSE)

do_pca(dds = dds_12, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"))
plot_vst_boxplot(transformed_12)

res_12 <- results_data(dds = dds_12, toCompare = "O", baseLevel = "Y", lfc = 0, alphaLim = 0.05, condition = "age", writeFile = FALSE)

res_12_1 <- results_data(dds = dds_12, toCompare = "O", baseLevel = "Y", lfc = 1, alphaLim = 0.05, condition = "age", writeFile = FALSE)

up_12 <- extract_DEG(res = res_12, direction = "up", lfc = 1, alphaLim = 0.05)
down_12 <- extract_DEG(res = res_12, direction = "down", lfc = 1, alphaLim = 0.05)


perform_GO_enrichment(gene_ids = up_12$ID, universe = rownames(counts_12))

extract_DEG_deseq2_data(res = res_12, direction = "up", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_up_12_deseq2_data")
extract_DEG_deseq2_data(res = res_12, direction = "down", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_down_12_deseq2_data")

extract_DEG_deseq2_data(res = res_12_1, direction = "up", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_up_12_deseq2_data_lfc_1")
extract_DEG_deseq2_data(res = res_12_1, direction = "down", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_down_12_deseq2_datalfc_1")


# 20-gene list mapped to Ensembl (GRCh38.p14)
gene_list <- data.frame(
  ID = c(
    "ENSG00000082781",  # ITGB5
    "ENSG00000117594",  # HSD11B1
    "ENSG00000112186",  # CAP2
    "ENSG00000162384",  # C1orf123
    "ENSG00000161057",  # PSMC2
    "ENSG00000115226",  # FNDC4
    "ENSG00000177905",  # COA3
    "ENSG00000126511",  # UBE2L6
    "ENSG00000153603",  # CYB5R3
    "ENSG00000160789",  # LMNA
    "ENSG00000169239",  # ZNF316
    "ENSG00000108344",  # CDK9
    "ENSG00000067083",  # LZTS2
    "ENSG00000046464",  # LTBR
    "ENSG00000141873",  # OCIAD2
    "ENSG00000077782",  # ITGA1
    "ENSG00000154016",  # MUM1L1
    "ENSG00000115354",  # CYC1
    "ENSG00000178279"   # GLIPR2
  ),
  stringsAsFactors = FALSE
)

do_MA_plot(res = res_12, gene_list = gene_list, fc = 2, fdr = 0.05)
