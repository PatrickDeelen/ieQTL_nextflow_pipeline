args <- commandArgs(trailingOnly = TRUE)

covars <- read.delim(args[1], sep = "\t", row.names = 1, as.is = T, check.names = F)
covar_name <- args[2]
output_prefix <- args[3]

pheno_is_factor <- F
if (length(unique(covars[,covar_name])) < 3) pheno_is_factor <- T

if (pheno_is_factor){
  for (level in unique(covars[,covar_name])){
    write.table(covars[covars[,covar_name] == level,], file = paste0(output_prefix, ".", level, ".txt"), sep = "\t", quote = F, col.names = NA)
  }
} else {
  q1 <- quantile(covars[,covar_name], probs = 0.25)
  q3 <- quantile(covars[,covar_name], probs = 0.75)
  
  write.table(covars[covars[,covar_name] < q1,], file = paste0(output_prefix, ".q1.txt"), sep = "\t", quote = F, col.names = NA)
  write.table(covars[covars[,covar_name] > q3,], file = paste0(output_prefix, ".q3.txt"), sep = "\t", quote = F, col.names = NA)
}
