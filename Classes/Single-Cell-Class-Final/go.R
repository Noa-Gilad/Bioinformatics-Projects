library(AnnotationDbi)
library(clusterProfiler)

over_expressed = read.delim("oocytes_gene_signatures_clusters.csv", header=TRUE, row.names = 1, sep = ",")

genes_to_test_gv <- over_expressed$GV_n[1:1000]

GO_results_BP_GV <- enrichGO(gene = genes_to_test_gv, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff=0.05)
as.data.frame(GO_results_BP_GV)
fit_up <- plot(barplot(GO_results_BP_GV, showCategory = 10))

genes_to_test_MII <- over_expressed$MII_n[1:1000]

GO_results_BP_MII <- enrichGO(gene = genes_to_test_MII, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff=0.05)
as.data.frame(GO_results_BP_MII)
fit_up <- plot(barplot(GO_results_BP_MII, showCategory = 10))
