###############################################################################
# Oocyte Aging – 75 Oocytes Dataset
# Author: Noa Gilad
#
# Description:
#   Differential expression and quality-control analysis of human oocyte
#   RNA-seq data from:
#
#     Llonch et al., Aging Cell (2021)
#
#   This dataset consists of 75 single human oocytes spanning young and
#   advanced maternal age groups. The analysis includes:
#     - Count preprocessing and QC
#     - PCA with and without batch correction (woman effect)
#     - Outlier detection and removal
#     - Differential expression using:
#         (1) GLMM (glmmSeq; age + random effect for woman)
#         (2) Averaged-per-woman DESeq2 with optional SVA correction
#
# Important Notes:
#   - This dataset was **explicitly excluded from merged / meta-analyses**
#     due to technical and biological differences relative to other aging
#     datasets.
#   - All analyses here are performed and interpreted independently.
#
# Inputs:
#   - count_75_FC.txt        : gene-level count matrix
#   - design_75.csv          : sample metadata (age, woman)
#
# Outputs:
#   - QC plots (histograms, boxplots, PCA)
#   - Differential expression result tables
#
# Dependencies:
#   - DESeq2
#   - glmmSeq
#   - Custom lab pipelines:
#       * deseq2_pipline.R
#       * preprocessing_and_QA.R
#       * glmm_pipeline.R
#
###############################################################################

# ~~~~~~~~~~~~ source pipeline ~~~~~~~~~~
source("C:\\Users\\gilad\\R\\pipelines\\deseq2_pipline.R", local = FALSE)
source("C:\\Users\\gilad\\R\\pipelines\\preprocessing_and_QA.R", local = FALSE)
source("C:\\Users\\gilad\\R\\pipelines\\glmm_pipeline.R", local = FALSE)

#  ~~~~~~~~~~~~~~~~~~~~~~~~~ analyse 75 oocytes paper data ~~~~~~~~~~~~~~~~~~~~~~~~~

df_75_oocytes <- read_tsv("count_75_FC.txt", comment = "#")
#df_75_clean <- df_75_oocytes %>% filter(!grepl("^__", gene_id))
df_75_clean <- df_75_oocytes %>% dplyr::select(Geneid, starts_with("SRR"))
colnames(df_75_clean) <- sub("Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(df_75_clean))
#df_75_clean_with_rownames <- df_75_clean %>% column_to_rownames(var = "gene_id")
df_75_clean_with_rownames <- df_75_clean %>% column_to_rownames(var = "Geneid")

colData_75_exp <- prep_coldata("design_75.csv", condition_col = "age", ref_level = "Y")

df_75_clean_with_rownames <- df_75_clean_with_rownames[, rownames(colData_75_exp)]
counts_75 <- prep_counts(counts = df_75_clean_with_rownames, readsLim = 15)
create_count_histogram(count_table = counts_75)

# ~~~~~~~~~~~~~~~~ QA ~~~~~~~~~~~~~~~~~~~
dds_75 <- DESeqDataSetFromMatrix(countData = counts_75, colData = colData_75_exp, design = ~ 1)
dds_75 <- estimateSizeFactors(dds_75)
normalized_counts_75 <- counts(dds_75, normalized = TRUE)
log_transformed_75 <- as.data.frame(log2(normalized_counts_75 + 1))
create_count_histogram(count_table = log_transformed_75, title = "log2 normalised counts per sample",  y_label = "log2 normalised counts")

boxplot(log_transformed_75, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,           # Make sample labels vertical
        col="lightblue", # Box color
        border="darkblue",
        main="RNA-Seq Expression Data",
        cex.axis=0.7,    # Reduce axis label size if needed
        par(mar=c(10,4,4,2) + 0.1)) # Increase bottom margin for labels

plotSampleViolin(log_transformed_75, sample_order = rownames(colData_75_exp))

do_pca(dds = dds_75, first_pc = 1, sec_pc = 2, intgroup = "age", ntop = 500, pca_colors = c("black","blue"), min_count = 100)
num_women_75 <- length(unique(colData(dds_75)$woman))
woman_colors_75 <- colorRampPalette(brewer.pal(8, "Dark2"))(num_women_75)
do_pca(dds = dds_75, first_pc = 1, sec_pc = 2, pca_colors = woman_colors_75, intgroup = c("woman"))

# Create a design matrix that preserves the age effect
design_mat_75 <- model.matrix(~ age, data = as.data.frame(colData(dds_75)))

# Extract the batch information from colData
woman_75_rm <- colData(dds_75)$woman

do_pca(dds = dds_75, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"), min_count = 10, batch = TRUE,  batch_rm = woman_75_rm, design = design_mat_75)

plotExpressionDistributions(log_transformed_75)

plot_all_pairs(log_counts = log_transformed_75, coldata = colData_75_exp, out_dir = "C:\\Users\\gilad\\R\\oocytes_aging\\75_results\\pairwise")


# ~~~~~~~~~~~~~~~~ outlier removal ~~~~~~~~~~~

# "SRR12746671", "SRR12746677", "SRR12746689", "SRR12746713" low counts 
outlier_samples_75 <- c("SRR12746671", "SRR12746677", "SRR12746689", "SRR12746713")

df_75_clean_with_rownames <- df_75_clean_with_rownames[, !(colnames(df_75_clean_with_rownames) %in% outlier_samples_75)]
counts_75 <- prep_counts(counts = df_75_clean_with_rownames, readsLim = 15)

colData_75_exp <- colData_75_exp[!(rownames(colData_75_exp) %in% outlier_samples_75), , drop = FALSE]


# ~~~~~~~~~~~~~~~~ Model fitting ~~~~~~~~~~~~~~~~~~

disp_75 <- estimate_disp(counts = counts_75, metadata = colData_75_exp)
size_75 <- estimate_size_factor(counts = counts_75)
formula_75 <- ~ age + ( 1 | woman)
fit_glmm_75 <- glmmSeq(modelFormula = formula_75, countdata = counts_75, metadata = colData_75_exp, dispersion = disp_75, sizeFactors = size_75, progress = TRUE)
fit_glmm_75 <- glmmQvals(fit_glmm_75)

glmm_res_75 <- formatGlmmSeqResults(object = fit_glmm_75, var = "age", target = "O")


up_75 <- extract_DE_genes(res = glmm_res_75, direction = "up", beta = 0, alphaLim = 0.05)
down_75 <- extract_DE_genes(res = glmm_res_75, direction = "down", beta = 0, alphaLim = 0.05)

extract_DE_data(res = glmm_res_75, direction = "up", beta = 0, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_up_75")
extract_DE_data(res = glmm_res_75, direction = "down", beta = 0, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_down_75")


# ~~~~~~~~~~~~~~~ second method - averaging and using deseq2 ~~~~~~~~~~~~~
params_75 <- sum_samples_by(coldata = colData_75_exp, count_table = counts_75, group_col = "woman", keep_cols = "age")
counts_ava_75 <- prep_counts(counts = params_75$count_table, readsLim = 15)
coldata_ava_75 <- params_75$coldata
coldata_ava_75[] <- lapply(coldata_ava_75, factor)
coldata_ava_75$age <- relevel(coldata_ava_75$age, ref = "Y")

dds_75_ava <- prepare_data(counts = counts_ava_75, coldata = coldata_ava_75, refer = "Y", des = ~ age, condition = "age")
normalized_counts_ave_75 <- counts(dds_75_ava, normalized = TRUE)
log_transformed_75_ave <- as.data.frame(log2(normalized_counts_ave_75 + 1))
create_count_histogram(count_table = log_transformed_75_ave, title = "log2 normalised counts per sample",  y_label = "log2 normalised counts")

transformed_lab_ava <-  vst(dds_75_ava, blind=FALSE)
do_pca(dds = dds_75_ava, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"), min_count = 50, min_samples = 2)

sva_out_75 <- apply_sva_correction(dds = dds_75_ava, full_model = ~ age, null_model = ~ 1, n_sv = 4, low_count_filter = 10)
ddssva_75   <- sva_out_75$ddssva
sv_vals_75  <- sva_out_75$sv_values 
mod_mat_75  <- model.matrix(~ age, colData(ddssva_75))
do_pca(
  ddssva_75,
  intgroup  = "age",
  sva       = TRUE,
  sva_rm    = sv_vals_75,
  design    = mod_mat_75
)


res_75_ava <- results_data(dds = ddssva_75, toCompare = "O", baseLevel = "Y", lfc = 0, alphaLim = 0.05, condition = "age", writeFile = FALSE)

up_75 <- extract_DEG(res = res_75_ava, direction = "up", lfc = 1, alphaLim = 0.05)
down_75 <- extract_DEG(res = res_75_ava, direction = "down", lfc = 1, alphaLim = 0.05)

extract_DEG_deseq2_data(res = res_75_ava, direction = "up", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_75_up_avg_sva")
extract_DEG_deseq2_data(res = res_75_ava, direction = "down", lfc = 1, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_75_up_down_avg_sva")
