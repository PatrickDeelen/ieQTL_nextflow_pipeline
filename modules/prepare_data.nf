#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


//gene_lengths = "$projectDir/data/Homo_sapiens.GRCh37.75.gene_lengths.txt.gz" 
//limix_annotation = "$projectDir/data/limix_gene_annotation_Ensembl71.txt.gz"

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
    path gene_lengths
    path limix_annotation
    val exp_platform

    output:
    path "*.TPM.txt.gz"

    script:
    """
    Rscript $projectDir/bin/normalize_TPM_and_rename.R \
        ${raw_expression} \
        ${gene_lengths} \
        $limix_annotation \
        expression.TPM.txt \
        $exp_platform
        
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
    val deconvolution_method
    val exptype

    output:
    path "cell_counts.txt"

    script:
    """
    Rscript $projectDir/bin/run_deconvolution.R \
        ${tpm_expression} \
        ${signature_matrix} \
	${deconvolution_method} \
        cell_counts.txt \
        ${exptype}

    """
}

/*
 * Combine major covariates with cell counts, genotype PCs and RNA-quality
 */
process CombineCovariatesRNAqual {
    publishDir params.outdir, mode: 'copy'
    tag "combine covars rnaqual"

    input:
    path general_covariates
    path cell_counts
    path genotype_PCs
    path gte
    path rna_qual

    output:
    path ("covariates.combined.txt"), emit: covariates_ch
    path ("*.distributions.pdf")

    script:
    """
    Rscript $projectDir/bin/combine_all_covariates.R -s ${general_covariates} -c ${cell_counts} -g ${genotype_PCs} -i ${gte} -o covariates.combined.txt -r ${rna_qual}
    ls -l ./
    """
}

/*
 * Combine major covariates with cell counts, genotype PCs
 */
process CombineCovariates {
    publishDir params.outdir, mode: 'copy'
    tag "combine covars rnaqual"

    input:
    path general_covariates
    path cell_counts
    path genotype_PCs
    path gte

    output:
    path "covariates.combined.txt"

    script:
    """
    Rscript $projectDir/bin/combine_all_covariates.R -s ${general_covariates} -c ${cell_counts} -g ${genotype_PCs} -i ${gte} -o covariates.combined.txt 
    """
}

/*
 * Creates a limix annotation file and a file with gene lengths from a GTF file
 */
process PrepareAnnotation {
    publishDir params.outdir, mode: 'copy'
    input:
    path gtf_annotation

    output:
    	path ("LimixAnnotationFile.txt"), emit: annotation_ch 
	path ("GeneLengths.txt"), emit: gene_lengths_ch
	path ("ChunkingFile.txt"), emit: chunks_ch

    script:
    """
	Rscript $projectDir/bin/createFeatureAnnotation.R --in_gtf ${gtf_annotation} --n_genes 500 --feature_name ENSG --out_dir ./
    """
}


process CalculateRNAQualityScore {
    input:
    path tmm_expression_data
    
    output:
	path ("RNA_quality.txt_CorrelationsWithAverageExpression.txt.gz"), emit: rnaquality_ch

    script:
    """
	python3 $projectDir/bin/correlate_samples_with_avg_gene_expression.py -ex ${tmm_expression_data} -op RNA_quality.txt -log2
    """

}

process NormalizeExpression {
    publishDir params.outdir, mode: 'copy'
    input:
      path(raw_expr)
      val(exp_platform)
      path(gte)
      path(check_sex)
      path(fam_file)

    output:
	path ('outputfolder_exp'), emit: expression_folder 
	path ('outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt'), emit: norm_expression_table
    
    shell:
    '''
       if [[ !{exp_platform} == "HT12v3" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_IlluminaHT12v3.txt
        elif [[ !{exp_platform} == "HuRef8" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_IlluminaHuRef8.txt
        elif [[ !{exp_platform} == "HT12v4" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_IlluminaHT12v4.txt
        elif [[ !{exp_platform} == "RNAseq" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_RNAseq.txt
        elif [[ !{exp_platform} == "AffyU219" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_AffyU219.txt
        elif [[ !{exp_platform} == "AffyHumanExon" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_AffyHumanExon.txt
        elif [[ !{exp_platform} == "RNAseq_HGNC" ]]; then
            probe_mapping_file=!{baseDir}/data/HgncToEnsemblProbeMatching.txt
        fi
        echo $probe_mapping_file
        
        outdir=${PWD}/outputfolder_exp/

        Rscript !{baseDir}/bin/ProcessExpression.R \
           -e !{raw_expr} \
           -l !{gte} \
           -p !{exp_platform} \
           -m $probe_mapping_file \
           -i !{check_sex} \
           -f !{fam_file} \
           -o ${outdir}    
    '''
}

process SplitCovariates {
    publishDir params.outdir, mode: 'copy'
    input:
    path normalized_expression_data
    path combined_covariates
    val covariate_name
    

    output:
        path "*txt"

    script:
    """
        Rscript $projectDir/bin/split_covariate_into_bins.R $combined_covariates $covariate_name $normalized_expression_data ./ 
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
        exp_type
	raw_expression_data
        normalized_expression_data
	signature_matrix_name
	deconvolution_method
        covariates
	gene_lengths
        limix_annotation
	genotype_pcs
	gte

    main:
	signature_matrix = "$projectDir/data/signature_matrices/" + signature_matrix_name  + ".txt.gz"
	if (exp_type == "RNAseq" || exp_type == "RNAseq_HGNC") {
	    cell_counts_ch = Deconvolution(TPM(raw_expression_data, gene_lengths, limix_annotation, exp_type), signature_matrix, deconvolution_method, exp_type)
            rnaquality_ch = CalculateRNAQualityScore(normalized_expression_data)
	    CombineCovariatesRNAqual(covariates,cell_counts_ch, genotype_pcs, gte, rnaquality_ch)
	    covariates_ch = CombineCovariatesRNAqual.out.covariates_ch
        } else {
	    cell_counts_ch = Deconvolution(normalized_expression_data, signature_matrix, deconvolution_method, exp_type)
	    covariates_ch = CombineCovariates(covariates,cell_counts_ch, genotype_pcs, gte)
	}
    emit:
        covariates_ch

}

