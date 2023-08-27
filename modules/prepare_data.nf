#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


gene_lengths = "$projectDir/data/Homo_sapiens.GRCh37.75.gene_lengths.txt.gz" 
limix_annotation = "$projectDir/data/limix_gene_annotation_Ensembl71.txt.gz"

/*
 * run TMM normalization on the raw expression file
 */
process TMM {
    tag "TMM normalization"
    //publishDir params.outdir, mode: 'copy'

    input:
    path raw_expression

    output:
    path "*TMM.txt.gz"

    script:
    """
    Rscript $projectDir/bin/normalizeTMM_cutoffCPM.R ${raw_expression} expression.TMM.txt 0.01
    gzip -f expression.TMM.txt
    """
}

/*
 * run TPM normalization on the raw expression file and rename gene ids to gene names
 */
process TPM {
    tag "TPM normalization"

    input:
    path raw_expression

    output:
    path "*.TPM.txt.gz"

    script:
    """
    Rscript $projectDir/bin/normalize_TPM_and_rename.R \
        ${raw_expression} \
        ${gene_lengths} \
        $limix_annotation \
        expression.TPM.txt
        
    gzip -f expression.TPM.txt
    """
}

/*
 * run deconvolution using nnls method with lm22 signature matrix
 */
process Deconvolution {
    tag "deconvolution"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path tpm_expression
    path signature_matrix

    output:
    path "cell_counts.txt"

    script:
    """
    Rscript $projectDir/bin/run_nnls_deconvolution.R \
        ${tpm_expression} \
        ${signature_matrix} \
        cell_counts.txt

    """
}

/*
 * run deconvolution using nnls method with LM22 signature matrix
 */
process CombineCovariates {
    tag "combine covars"
    
    input:
    path cell_counts
    path covariates

    output:
    path "covariates.INT.txt"

    script:
    """
    Rscript $projectDir/bin/combine_covariates_run_INT.R ${covariates} ${cell_counts} covariates.INT.txt
    
    """
}

workflow TMM_TRANSFORM_EXPRESSION {
    take:
        raw_expression_data

    main:
        tmm_expression_ch = TMM(raw_expression_data)
        
    emit:
        tmm_expression_ch

}

workflow PREPARE_COVARIATES {
    take:
        raw_expression_data
        signature_matrix
        covariates

    main:
        cell_counts_ch = Deconvolution(TPM(raw_expression_data), signature_matrix)
        covariates_ch = CombineCovariates(cell_counts_ch, covariates)
        
    emit:
        covariates_ch

}


