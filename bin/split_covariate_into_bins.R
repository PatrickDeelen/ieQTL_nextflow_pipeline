library(data.table)

args <- commandArgs(trailingOnly = TRUE)

covars <- read.delim(args[1], sep = "\t", row.names = 1, as.is = T, check.names = F)
covar_name <- args[2]
norm_expression <- read.delim(args[3], sep = "\t", row.names = 1, as.is = T, check.names = F)
output_prefix <- args[4]

pheno_is_factor <- F
if (length(unique(covars[,covar_name])) < 3) pheno_is_factor <- T


exp_summary <- function(x){
  per_gene_mean <- apply(x, 1, mean)
  per_gene_median <- apply(x, 1, median)
  per_gene_min <- apply(x, 1, min)
  per_gene_max <- apply(x, 1, max)
  per_gene_sd <- apply(x, 1, sd)
  nr_values <- apply(x, 1, function(x) length(x))
  unique_values <- apply(x, 1, function(x) length(unique(x)))
  zero_values <- apply(x, 1, function(x) length(x[x == 0]))
  per_gene_shapiro <- apply(x, 1, shap_test)
  
  gene_summary <- data.table(gene = rownames(x), 
                             mean = per_gene_mean,
                             median = per_gene_median,
                             min = per_gene_min,
                             max = per_gene_max,
                             sd = per_gene_sd,
                             nr_values = nr_values,
                             nr_unique_values = unique_values,
                             nr_of_zero_values = zero_values,
                             shapiro_P = per_gene_shapiro)
  
  return(gene_summary)
}

shap_test <- function(x){
  if (length(unique(x)) > 1) {
    return(shapiro.test(x)$p.value)
  } else {
    return(NA)
  }
}

if (pheno_is_factor){
  for (level in unique(covars[,covar_name])){
    covar_for_bin <- covars[covars[,covar_name] == level,]
    expr_summary_for_bin <- exp_summary(norm_expression[row.names(covar_for_bin), ])
    
    write.table(covar_for_bin, file = paste0(output_prefix, "covariates.", covar_name, "_", level, ".txt"), sep = "\t", quote = F, col.names = NA)
    write.table(expr_summary_for_bin, file = paste0(output_prefix, "expression_summary.", covar_name, "_", level, ".txt"), sep = "\t", quote = F, col.names = NA)
    
  }
} else {
  q1 <- quantile(covars[,covar_name], probs = 0.25)
  q3 <- quantile(covars[,covar_name], probs = 0.75)
  covars_q1 <- covars[covars[,covar_name] < q1,]
  covars_q3 <- covars[covars[,covar_name] > q3,]
  
  expr_summary_q1 <- exp_summary(norm_expression[row.names(covars_q1), ])
  expr_summary_q3 <- exp_summary(norm_expression[row.names(covars_q3), ])
  
  write.table(covars_q1, file = paste0(output_prefix, "covariates.", covar_name, ".q1.txt"), sep = "\t", quote = F, col.names = NA)
  write.table(covars_q3, file = paste0(output_prefix, "covariates.", covar_name, ".q3.txt"), sep = "\t", quote = F, col.names = NA)

  write.table(expr_summary_q1, file = paste0(output_prefix, "expression_summary.", covar_name, ".q1.txt"), sep = "\t", quote = F, col.names = NA)
  write.table(expr_summary_q3, file = paste0(output_prefix, "expression_summary.", covar_name, ".q3.txt"), sep = "\t", quote = F, col.names = NA)
}
                                                                                                                                                                                                                   
