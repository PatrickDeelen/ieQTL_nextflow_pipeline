args <- commandArgs(trailingOnly = TRUE)

#setwd("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/data")
bulk_expr <- as.matrix(read.delim(args[1], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
gene_lengths <- read.delim(args[2], row.names = 1, header = F, as.is = T, check.names = F, sep = "\t")
probe_annotation <- read.delim(args[3], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t")
outfile=args[4]

#bulk_expr <- as.matrix(read.delim(gzfile("gene.counts-LLDBIOSSamples.LLD_subset.txt.gz"), row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
#gene_lengths <- read.delim("Homo_sapiens.GRCh37.75.gene_lengths.txt.gz", row.names = 1, header = F, as.is = T, check.names = F, sep = "\t")
#probe_annotation <- read.delim("ProbeAnnotation_STARv2.3.0e_Ensembl71.txt.gz", row.names = 2, header = T, as.is = T, check.names = F, sep = "\t")

# TPM normalize
bulk_expr <- bulk_expr[which(rowSums(bulk_expr) != 0),]
shared_genes <- intersect(row.names(bulk_expr), row.names(gene_lengths))
cat("Number of genes in expr matrix: ", nrow(bulk_expr), "\nNumber of genes shared with gene length table: ", length(shared_genes), "\n")
gene_lengths <- gene_lengths[shared_genes,]
bulk_expr <- bulk_expr[shared_genes,]

#bulk_expr_tpm <- get_TPM(bulk_expr[shared_genes,], gene_lengths)
bulk_expr_tpm <- do.call(cbind, lapply(1:ncol(bulk_expr), function(i) {
  rate = log(bulk_expr[,i]) - log(gene_lengths)
  denom = log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}))
colnames(bulk_expr_tpm) <- colnames(bulk_expr)

# Convert Ensembl ids to gene names
new_names <-  probe_annotation[match(row.names(bulk_expr_tpm), row.names(probe_annotation), nomatch = "NA"),"GeneName"]
new_names <- gsub(" ", "", new_names, fixed = TRUE)
row.names(bulk_expr_tpm) <- new_names

n_occur <- data.frame(table(new_names))
dups <- n_occur[n_occur$Freq > 1,"new_names"]
bulk_expr_tpm <- bulk_expr_tpm[!row.names(bulk_expr_tpm) %in% dups,]

write.table(na.omit(bulk_expr_tpm), file = args[4], sep = "\t", quote = F, col.names = NA)
#write.table(bulk_expr_tpm, file = paste0("gene.counts-LLDBIOSSamples.LLD_subset.txt.gz", ".TPM.txt"), sep = "\t", quote = F, col.names = NA)



