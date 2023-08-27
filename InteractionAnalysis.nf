#! /usr/bin/env nextflow


nextflow.enable.dsl = 2


gene_lengths = "$projectDir/data/Homo_sapiens.GRCh37.75.gene_lengths.txt.gz" 
limix_annotation = "$projectDir/data/limix_gene_annotation_Ensembl71.txt.gz"
signature_matrix = "$projectDir/data/signature_matrices/LM22.txt.gz"
outdir = params.outdir

raw_expr_ch = Channel.fromPath(params.expfile)
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)

include { TMM_TRANSFORM_EXPRESSION; PREPARE_COVARIATES } from './modules/prepare_data.nf'

workflow {
    tmm_ch = TMM_TRANSFORM_EXPRESSION(raw_expr_ch)
    covariates_ch = PREPARE_COVARIATES(raw_expr_ch, signature_matrix, covars_ch)

    //cell_counts_ch.flatten().collectFile(name: 'cell_counts.txt', keepHeader: true, sort: true, storeDir: "${params.outdir}")
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
