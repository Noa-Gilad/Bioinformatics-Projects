###############################################################################
# Oocyte Aging – Lab Data Analysis
# Author: Noa Gilad
#
# Description:
#   Differential expression analysis of human oocyte aging using in-house
#   lab-generated RNA-seq data. This script performs quality assessment,
#   outlier detection, and differential expression using two complementary
#   approaches:
#
#     1. GLMM-based modeling (glmmSeq) accounting for repeated measures
#        from the same woman.
#     2. Aggregation by donor followed by DESeq2 analysis, with optional
#        SVA correction.
#
#   The analysis compares old (O) versus young (Y) oocytes and evaluates
#   robustness across modeling strategies.
#
# Input:
#   - Count table generated from HTSeq/featureCounts (e.g. count_lab_FC.txt)
#   - Experimental design file (design_lab.csv)
#
# Outputs:
#   - QC plots (histograms, boxplots, PCA)
#   - Differential expression result tables for both methods
#
# Dependencies:
#   - Custom pipelines:
#       * glmm_pipeline.R
#       * deseq2_pipline.R
#       * preprocessing_and_QA.R
#   - Packages: DESeq2, glmmSeq, sva, tidyverse, ggplot2
#
# Notes:
#   - Paths are user-specific and should be adapted locally.
#   - This script is intended for reproducibility and documentation,
#     not as a standalone package.
###############################################################################

# ~~~~~~~~~~~~ source pipeline ~~~~~~~~~~
source("C:\\Users\\gilad\\R\\pipelines\\glmm_pipeline.R", local = FALSE)
source("C:\\Users\\gilad\\R\\pipelines\\deseq2_pipline.R", local = FALSE)
source("C:\\Users\\gilad\\R\\pipelines\\preprocessing_and_QA.R", local = FALSE)

#  ~~~~~~~~~~~~~~~~~~~~~~~~~ analyse lab data ~~~~~~~~~~~~~~~~~~~~~~~~~

df_lab <- read_tsv("count_lab_FC.txt", comment = "#")
# df_lab_clean <- df_lab %>% filter(!grepl("^__", gene_id))
df_lab_clean <- df_lab %>% dplyr::select(Geneid, starts_with("R_"))
colnames(df_lab_clean) <- sub("_dedup\\.bam$", "", colnames(df_lab_clean))
df_lab_clean_with_rownames <- df_lab_clean %>% column_to_rownames(var = "Geneid") # for htseq - gene_id


colData_lab_exp <- prep_coldata("design_lab.csv", condition_col = "age", ref_level = "Y")

df_lab_clean_with_rownames <- df_lab_clean_with_rownames[, rownames(colData_lab_exp)]
counts_lab <- prep_counts(counts = df_lab_clean_with_rownames, readsLim = 20)
create_count_histogram(count_table = counts_lab)

# ~~~~~~~~~~~~~~~~ QA ~~~~~~~~~~~~~~~~~~~
dds_lab <- DESeqDataSetFromMatrix(countData = counts_lab, colData = colData_lab_exp, design = ~ 1)
dds_lab <- estimateSizeFactors(dds_lab)
normalized_counts <- counts(dds_lab, normalized = TRUE)
log_transformed_lab <- as.data.frame(log2(normalized_counts + 1))
create_count_histogram(count_table = log_transformed_lab, title = "log2 normalised counts per sample",  y_label = "log2 normalised counts")

boxplot(log_transformed_lab, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,           # Make sample labels vertical
        col="lightblue", # Box color
        border="darkblue",
        main="RNA-Seq Expression Data",
        cex.axis=0.7,    # Reduce axis label size if needed
        par(mar=c(10,4,4,2) + 0.1)) # Increase bottom margin for labels

plotSampleViolin(log_transformed_lab, sample_order = rownames(colData_lab_exp))

do_pca(dds = dds_lab,first_pc = 1, sec_pc = 2, intgroup = "age", ntop = 500, pca_colors = c("black","blue"))
num_women <- length(unique(colData(dds_lab)$woman))
woman_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(num_women)
do_pca(dds = dds_lab, first_pc = 1, sec_pc = 2, pca_colors = woman_colors, intgroup = c("woman"))


# Create a design matrix that preserves the age effect
design_lab <- model.matrix(~ age, data = as.data.frame(colData(dds_lab)))

# Extract the batch information from colData
woman <- colData(dds_lab)$woman
 
do_pca(dds = dds_lab,first_pc = 1, sec_pc = 2, intgroup = "age", ntop = 500, pca_colors = c("black","blue"), min_count = 15 , batch = TRUE,  batch_rm = woman, design = design_lab)

plot_all_pairs(log_counts = log_transformed_lab, coldata = colData_lab_exp, out_dir = "C:\\Users\\gilad\\R\\oocytes_aging\\lab_results\\pairwise")


# ~~~~~~~~~~~~~~~~ outlier removal ~~~~~~~~~~~

outlier_samples_lab <- c("R_2b_S36")
df_lab_clean_with_rownames <- df_lab_clean_with_rownames[, !(colnames(df_lab_clean_with_rownames) %in% outlier_samples_lab)]
counts_lab <- prep_counts(counts = df_lab_clean_with_rownames, readsLim = 15)
colData_lab_exp <- colData_lab_exp[!(rownames(colData_lab_exp) %in% outlier_samples_lab), , drop = FALSE]


# ~~~~~~~~~~~~~~~~ Model fitting ~~~~~~~~~~~~~~~~~~
disp <- estimate_disp(counts = counts_lab, metadata = colData_lab_exp)
size <- estimate_size_factor(counts = counts_lab)
formula <- ~ age + ( 1 | woman)
fit_glmm <- glmmSeq(modelFormula = formula, countdata = counts_lab, metadata = colData_lab_exp, dispersion = disp, sizeFactors = size, progress = TRUE)
fit_glmm <- glmmQvals(fit_glmm)

glmm_res <- formatGlmmSeqResults(object = fit_glmm, var = "age", target = "O")


up_lab <- extract_DE_genes(res = glmm_res, direction = "up", beta = 0, alphaLim = 0.05)
down_lab <- extract_DE_genes(res = glmm_res, direction = "down", beta = 0, alphaLim = 0.05)

extract_DE_data(res = glmm_res, direction = "up", beta = 0, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_up_lab_data")
extract_DE_data(res = glmm_res, direction = "down", beta = 0, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_down_lab_data")



# ~~~~~~~~~~~~~~~ second method - averaging and using deseq2 and SVA ~~~~~~~~~~~~~
params <- sum_samples_by(coldata = colData_lab_exp, count_table = counts_lab, group_col = "woman", keep_cols = "age")  # either averaging or summing doesn't result in significant genes
counts_ava_lab <- prep_counts(counts = params$count_table, readsLim = 15)
coldata_ava_lab <- params$coldata
coldata_ava_lab[] <- lapply(coldata_ava_lab, factor)
coldata_ava_lab$age <- relevel(coldata_ava_lab$age, ref = "Y")

dds_lab_ava <- prepare_data(counts = counts_ava_lab, coldata = coldata_ava_lab, refer = "Y", des = ~ age, condition = "age")
normalized_counts_ave <- counts(dds_lab_ava, normalized = TRUE)
log_transformed_lab_ave <- as.data.frame(log2(normalized_counts_ave + 1))
create_count_histogram(count_table = log_transformed_lab_ave, title = "log2 normalised counts per sample",  y_label = "log2 normalised counts")

transformed_lab_ava <-  vst(dds_lab_ava, blind=FALSE)
do_pca(dds = dds_lab_ava, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"), min_count = 50, min_samples = 2)

sva_out <- apply_sva_correction(dds = dds_lab_ava, full_model = ~ age, null_model = ~ 1, n_sv = 4, low_count_filter = 10)
ddssva   <- sva_out$ddssva
sv_vals  <- sva_out$sv_values 
mod_mat  <- model.matrix(~ age, colData(ddssva))
do_pca(
  ddssva,
  intgroup  = "age",
  sva       = TRUE,
  sva_rm    = sv_vals,
  design    = mod_mat
)


res_lab_ava <- results_data(dds = ddssva, toCompare = "O", baseLevel = "Y", lfc = 0, alphaLim = 0.05, condition = "age", writeFile = FALSE)

up_lab <- extract_DEG(res = res_lab_ava, direction = "up", lfc = 1, alphaLim = 0.05)
down_lab <- extract_DEG(res = res_lab_ava, direction = "down", lfc = 1, alphaLim = 0.05)

extract_DEG_deseq2_data(res = res_lab_ava, direction = "up", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_lab_up_avg_sva")
extract_DEG_deseq2_data(res = res_lab_ava, direction = "down", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_lab_down_avg_sva")
