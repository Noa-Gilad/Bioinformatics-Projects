# ==============================================================================
# Oocyte aging – merged analysis (lab data + 12-oocyte dataset)
# Author: Noa Gilad
#
# Description:
#   Differential expression and quality-control analysis of a merged human
#   oocyte aging dataset combining:
#     1. In-house lab oocyte RNA-seq data
#     2. Public 12-oocyte dataset (Yuan et al.)
#
#   The analysis evaluates age-associated transcriptional changes while
#   explicitly accounting for dataset origin and donor (woman) effects.
#
# Key analyses:
#   - Count matrix merging by Geneid
#   - QC and exploratory analysis (normalization, PCA, batch inspection)
#   - GLMM-based differential expression:
#       age + dataset + (1 | woman)
#   - Alternative strategy using per-woman averaging + DESeq2
#   - Batch correction using explicit covariates and SVA
#
# Important notes:
#   - The 75-oocyte dataset (Llonch et al., Aging Cell 2021) is intentionally
#     excluded from this merged analysis due to incompatibility identified
#     during cross-dataset QC.
#   - Dataset identity is treated as a biological/technical covariate and
#     explicitly modeled where applicable.
#
# References:
#   - Yuan et al., human oocyte aging (12 oocytes)
#   - In-house lab oocyte aging dataset
#
# ==============================================================================

# ~~~~~~~~~~~~ source pipeline ~~~~~~~~~~
source("C:\\Users\\gilad\\R\\pipelines\\deseq2_pipline.R", local = FALSE)
source("C:\\Users\\gilad\\R\\pipelines\\preprocessing_and_QA.R", local = FALSE)
source("C:\\Users\\gilad\\R\\pipelines\\glmm_pipeline.R", local = FALSE)

#  ~~~~~~~~~~~~~~~~~~~~~~~~~ analyse merged data (lab data and the 12 oocytes dataset) ~~~~~~~~~~~~~~~~~~~~~~~~~

# Merge the data frames by the common "GeneID" column using full_join
df_merged <- df_lab_clean %>%
  full_join(df_12_clean, by = "Geneid") # change to gene_id for htseq

# Write the merged data to a TSV file
write_csv(df_merged, "merged_counts_12_lab.csv")

df_merged <- df_merged %>% column_to_rownames(var = "Geneid")

colData_merged <- prep_coldata("merged_lab12_design.csv", condition_col = "age", ref_level = "Y")

colData_merged <- colData_merged[order(colData_merged$age, decreasing = TRUE),]
df_merged <- df_merged[, rownames(colData_merged)]

counts_merged <- prep_counts(counts = df_merged, readsLim = 33)

dds_merged <- DESeqDataSetFromMatrix(countData = counts_merged, colData = colData_merged, design = ~ 1)
dds_merged <- estimateSizeFactors(dds_merged)
normalized_counts_merged <- counts(dds_merged, normalized = TRUE)
log_transformed_merged <- as.data.frame(log2(normalized_counts_merged + 1))
create_count_histogram(count_table = log_transformed_merged, title = "log2 normalised counts per sample",  y_label = "log2 normalised counts")

boxplot(log_transformed_lab, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,           # Make sample labels vertical
        col="lightblue", # Box color
        border="darkblue",
        main="RNA-Seq Expression Data",
        cex.axis=0.7,    # Reduce axis label size if needed
        par(mar=c(10,4,4,2) + 0.1)) # Increase bottom margin for labels


do_pca(dds = dds_merged ,first_pc = 1, sec_pc = 2, intgroup = "age", ntop = 500, pca_colors = c("black","blue"))
do_pca(dds = dds_merged ,first_pc = 1, sec_pc = 2, intgroup = "dataset", ntop = 500, pca_colors = c("black", "red"))
num_women_merged <- length(unique(colData(dds_merged)$woman))
woman_colors_merged <- colorRampPalette(brewer.pal(8, "Dark2"))(num_women_merged)
do_pca(dds = dds_merged, first_pc = 1, sec_pc = 2, pca_colors = woman_colors_merged, intgroup = c("woman"))

# Create a design matrix that preserves the age effect
design_mat <- model.matrix(~ age, data = as.data.frame(colData(dds_merged)))

# Extract the batch information from colData
batch <- colData(dds_merged)$dataset

do_pca(dds = dds_merged ,first_pc = 1, sec_pc = 2, intgroup = "age", ntop = 500, pca_colors = c("black","blue"), min_count = 15 , batch = TRUE,  batch_rm = batch, design = design_mat)

outlier_samples_m <- c("SRR12330964", "SRR12330963", "SRR12330965", "R_2b_S36")
df_merged <- df_merged[, !(colnames(df_merged) %in% outlier_samples_m)]
counts_merged <- prep_counts(counts = df_merged, readsLim = 30)
colData_merged <- colData_merged[!(rownames(colData_merged) %in% outlier_samples_m), , drop = FALSE]

disp_m <- estimate_disp(counts = counts_merged, metadata = colData_merged)
size_m <- estimate_size_factor(counts = counts_merged)
formula_merged <- ~ age + dataset + (1 | woman)
fit_glmm_m <- glmmSeq(modelFormula = formula_merged, countdata = counts_merged, metadata = colData_merged, dispersion = disp_m, sizeFactors = size_m, progress = TRUE)
fit_glmm_m <- glmmQvals(fit_glmm_m)

glmm_res_m <- formatGlmmSeqResults(object = fit_glmm_m, var = "age", target = "O")


up_merged <- extract_DE_genes(res = glmm_res_m, direction = "up", beta = 0, alphaLim = 0.05)
down_merged <- extract_DE_genes(res = glmm_res_m, direction = "down", beta = 0, alphaLim = 0.05)

extract_DE_data(res = glmm_res_m, direction = "up", beta = 0, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_up_12_and_lab")
extract_DE_data(res = glmm_res_m, direction = "down", beta = 0, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_down_12_and_lab")

# ~~~~~~~~~~~~~~~ second method - averaging and using deseq2 ~~~~~~~~~~~~~
params_m <- sum_samples_by(coldata = colData_merged, count_table = counts_merged, group_col = "woman", keep_cols = c("age", "dataset"))
counts_ave_merged <- prep_counts(counts = params_m$count_table, readsLim = 25)
coldata_ave_merged <- params_m$coldata
coldata_ave_merged[] <- lapply(coldata_ave_merged, factor)
coldata_ave_merged$age <- relevel(coldata_ave_merged$age, ref = "Y")

dds_merged_ave <- prepare_data(counts = counts_ave_merged, coldata = coldata_ave_merged, refer = "Y", des = ~ age + dataset, condition = "age")
normalized_counts_ave_merged <- counts(dds_merged_ave, normalized = TRUE)
log_transformed_merged_ave <- as.data.frame(log2(normalized_counts_ave_merged + 1))
create_count_histogram(count_table = log_transformed_merged_ave, title = "log2 normalised counts per sample",  y_label = "log2 normalised counts")

transformed_merged_ave <-  vst(dds_merged_ave, blind=FALSE)
do_pca(dds = dds_merged_ave, first_pc = 1, sec_pc = 2, pca_colors = c("black","blue"), intgroup =  c("age"), min_count = 50, min_samples = 2, batch = TRUE,
       batch_rm = colData(dds_merged_ave)$dataset,
       design    = model.matrix(~ age, colData(dds_merged_ave)))

sva_out_merged <- apply_sva_correction(dds = dds_merged_ave, full_model = ~ age + dataset, null_model = ~ dataset, n_sv = 2, low_count_filter = 10)
ddssva_merged   <- sva_out_merged$ddssva
sv_vals_merged  <- sva_out_merged$sv_values 
mod_mat_merged  <- model.matrix(~ age, colData(ddssva_merged))
do_pca(
  ddssva_merged,
  intgroup  = "age",
  sva       = TRUE,
  sva_rm    = sv_vals_merged,
  batch = TRUE,
  batch_rm = colData(ddssva_merged)$dataset,
  design    = mod_mat_merged
)

res_merged_ave <- results_data(dds = ddssva_merged, toCompare = "O", baseLevel = "Y", lfc = 2, alphaLim = 0.05, condition = "age", writeFile = FALSE)

up_merged_ave <- extract_DEG(res = res_merged_ave, direction = "up", lfc = 2, alphaLim = 0.05)
down_merged_ave <- extract_DEG(res = res_merged_ave, direction = "down", lfc = 2, alphaLim = 0.05)

extract_DEG_deseq2_data(res = res_merged_ave, direction = "up", lfc = 2, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_12lab_up_avg_sva2")
extract_DEG_deseq2_data(res = res_merged_ave, direction = "down", lfc = 2, alphaLim = 0.05, writeFile = TRUE, title = "O_vs_Y_12lab_down_avg_sva2")
