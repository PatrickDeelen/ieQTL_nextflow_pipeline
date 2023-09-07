#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

//signature_matrix = "$projectDir/data/signature_matrices/LM22.txt.gz"

params.qtls_to_test = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input//qtls_to_test.txt"
params.bfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input//LLD_genotypes_flt"
params.gte = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/LLD_gte.txt"
params.outdir = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/"
params.covariate_to_test = "age"
params.num_perm = 0
params.gtf = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/Homo_sapiens.GRCh37.75.gtf.gz"

params.signature_matrix_name = "LM22"
params.deconvolution_method = "nnls"


raw_expr_ch = Channel.fromPath(params.expfile)
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)
gtf_annotation_ch = Channel.fromPath(params.gtf)

Channel
    .from(params.bfile)
    .ifEmpty { exit 1, "Input plink prefix not found!" }
    .map { genotypes -> [file("${genotypes}.bed"), file("${genotypes}.bim"), file("${genotypes}.fam")]}
    .set { bfile_ch }


include { TMM_TRANSFORM_EXPRESSION; PREPARE_COVARIATES; prepare_annotation } from './modules/prepare_data.nf'
include { ieQTL_mapping } from './modules/interaction_analysis.nf'

workflow {
    tmm_ch = TMM_TRANSFORM_EXPRESSION(raw_expr_ch)

    prepare_annotation(gtf_annotation_ch)

    prepare_annotation.out.chunks_ch.splitCsv( header: false ).map { row -> row[0] }.set { chunk_ch }
    covariates_ch = PREPARE_COVARIATES(raw_expr_ch, params.signature_matrix_name, params.deconvolution_method,covars_ch, prepare_annotation.out.gene_lengths_ch, prepare_annotation.out.annotation_ch)

    eqtl_ch = tmm_ch.combine(bfile_ch).combine(covariates_ch).combine(prepare_annotation.out.annotation_ch).combine(Channel.fromPath(params.gte)).combine(Channel.fromPath(params.qtls_to_test)).combine(Channel.of(params.covariate_to_test)).combine(chunk_ch)
    
    results_ch = ieQTL_mapping(eqtl_ch)

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
