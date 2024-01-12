#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="DataQc"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load PLINK
ml nextflow

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR=../singularitycache
export NXF_HOME=../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

# Define paths
c=LLS_OmniExpr
# Genotype data
#[full path to the folder with imputed filtered vcf files produced by eQTLGen pipeline 2_Imputation step (postimpute folder)]
vcf_dir_path=/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/2_Imputation/out/${c}/postimpute/
raw_exp_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/data/${c}/${c}_raw_expression.txt.gz
gte_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/data/${c}/${c}.gte
norm_exp_path=/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/1_DataQC/out/${c}/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt
covariate_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/BIOS_covariates.txt
genotype_pcs_path=/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/1_DataQC/out/${c}/outputfolder_gen/gen_PCs/GenotypePCs.txt

exp_platform=RNAseq
cohort_name=$c
genome_build="GRCh38"
covariate_to_test=gender_F1M2
genes_to_test=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/bios_sign_genes.txt
output_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/results/${c}/

# Additional settings and optional arguments for the command


# Command:
NXF_VER=21.10.6 nextflow run /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/ieQTL_nextflow_pipeline/InteractionAnalysis.nf \
--vcf_dir ${vcf_dir_path} \
--raw_expfile ${raw_exp_path} \
--norm_expfile ${norm_exp_path} \
--gte ${gte_path} \
--covariates $covariate_path \
--exp_platform ${exp_platform} \
--cohort_name ${cohort_name} \
--covariate_to_test $covariate_to_test \
--genes_to_test $genes_to_test \
--genotype_pcs $genotype_pcs_path \
--outdir ${output_path}  \
-profile slurm,singularity \
-resume