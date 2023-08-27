args <- commandArgs(trailingOnly = TRUE)
library(nnls)

model_nnls = function(y, x){
  fit <- apply(y, 2, function(z) {
    nnls(A = x, b = as.vector(z))})
  
  coeff <- data.frame(lapply(fit, function(z) z$x))
  colnames(coeff) <- colnames(y)
  rownames(coeff) <- colnames(x)
  coeff <- t(coeff)
  
  return(coeff)
}

expr_fname = args[1]
sign_fname = args[2]
out_fname = args[3]

bulk_expr_tpm <- as.matrix(read.delim(expr_fname, row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
lm22 <- as.matrix(read.delim(sign_fname, row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))

shared_genes <- intersect(row.names(bulk_expr_tpm), row.names(lm22))
cat(length(shared_genes), "genes shared between expression table and signature matrix\n")

coef <- model_nnls(x = lm22[shared_genes,], y = bulk_expr_tpm[shared_genes,])
proportions <- round((100/max(rowSums(coef), na.rm=TRUE))*coef, 4)

write.table(proportions, file = out_fname, sep = "\t", quote = F, col.names = NA)

