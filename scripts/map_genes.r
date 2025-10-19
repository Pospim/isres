library(biomaRt)

mart <- useMart("ensembl", dataset = "ggallus_gene_ensembl")

gene_df <- read.csv("~/Desktop/bca/chicken_isgs_mapped.csv",
                     stringsAsFactors = FALSE,
                     sep = "\t")
gene_ids <- gene_df$ensembl_gene_id

mapping <- getBM(attributes = c("ensembl_gene_id", "refseq_mrna"),
                 filters = "ensembl_gene_id",
                 values = gene_ids, mart = mart)

head(mapping)

write.table(mapping, file = "~/Desktop/bca/chicken_isgs_mapped2.csv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
