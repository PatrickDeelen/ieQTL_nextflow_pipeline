#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


gene_lengths = "$projectDir/data/Homo_sapiens.GRCh37.75.gene_lengths.txt.gz" 
limix_annotation = "$projectDir/data/limix_gene_annotation_Ensembl71.txt.gz"
signature_matrix = "$projectDir/data/signature_matrices/LM22.txt.gz"
//outdir = params.outdir

baseDir="/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/"
params.qtls_to_test = "${baseDir}/qtls_to_test.txt"
params.bfile = "${baseDir}/LLD_genotypes_flt"
params.gte = "${baseDir}/LLD_gte.txt"
params.outdir = "${baseDir}/../output/"
params.covariate_to_test = "gender_F1M2"
params.num_perm = 0

//Q: do I need to explicitly initialize channels for each input parameter?
raw_expr_ch = Channel.fromPath(params.expfile)
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)

include { TMM_TRANSFORM_EXPRESSION; PREPARE_COVARIATES } from './modules/prepare_data.nf'
include { ieQTL_mapping } from './modules/interaction_analysis.nf'

workflow {
    tmm_ch = TMM_TRANSFORM_EXPRESSION(raw_expr_ch)
    covariates_ch = PREPARE_COVARIATES(raw_expr_ch, signature_matrix, covars_ch)
    results_ch = ieQTL_mapping(tmm_ch, covariates_ch, limix_annotation, params.outdir, params.covariate_to_test, params.num_perm)
    
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
