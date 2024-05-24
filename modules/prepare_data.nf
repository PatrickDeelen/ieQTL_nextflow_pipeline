#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * run TPM normalization on the raw expression file and rename gene ids to gene names
 */
process TPM {
    tag "TPM normalization"
    label "short"

    input:
    path raw_expression
    path gene_lengths


    output:
    path "*.TPM.txt.gz"

    script:
    """
    Rscript $projectDir/bin/normalize_TPM.R \
        ${raw_expression} \
        ${gene_lengths} \
        expression.TPM.txt 
        
    gzip -f expression.TPM.txt
    """
}

/*
 * run deconvolution using nnls method with lm22 signature matrix
 */
process Deconvolution {
    tag "deconvolution"
    label "medium1"

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


process TransformCovariates {
    echo true
    label "medium1"

    publishDir params.outdir, mode: 'copy'

    input:
    val cov_type
    path covariate_file
    path gte
    val rename
    val INT
    val cut

    output:
    path "*.transformed.txt"
    

    script:
    """
        Rscript $projectDir/bin/transform_covariates.R \
        -i ${covariate_file} \
        -g ${gte} \
        --rename ${rename} \
        --int ${INT} \
        -o ${cov_type}.transformed.txt \
        --cut ${cut}

        cp *pdf ${params.outdir}/
    """
}


/*
 * Combine major covariates with cell counts, genotype PCs and RNA-quality
 */
process CombineCovariatesRNAqual {
    label "medium1"

    publishDir params.outdir, mode: 'copy'

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
    if (cell_counts != "NA")
      """
      Rscript $projectDir/bin/combine_all_covariates.R -s ${general_covariates} -c ${cell_counts} -g ${genotype_PCs} -i ${gte} -o covariates.combined.txt -r ${rna_qual}
      """
    else
      """
      Rscript $projectDir/bin/combine_all_covariates.R -s ${general_covariates} -g ${genotype_PCs} -i ${gte} -o covariates.combined.txt -r ${rna_qual}
      """
}

/*
 * Combine major covariates with cell counts, genotype PCs
 */
process CombineCovariates {
    label "medium1"

    publishDir params.outdir, mode: 'copy'

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

process CalculateRNAQualityScore {
    label "short"
    
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
    label "medium2"
    echo true
    publishDir params.outdir, mode: 'copy'
    input:
      path(raw_expr)
      path(norm_expr)
      val(exp_platform)
      path(gte)

    output:
	path ('outputfolder_exp'), emit: expression_folder 
	path ('outputfolder_exp/exp_data_preprocessed.txt'), emit: norm_expression_table
    
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
                
        outdir=${PWD}/outputfolder_exp/

        Rscript !{baseDir}/bin/ProcessExpression_v2.R \
           -e !{raw_expr} \
           -n !{norm_expr} \
           -l !{gte} \
           -p !{exp_platform} \
           -m $probe_mapping_file \
           -o ${outdir}    
    '''
}
process SplitCovariates {
    label "medium1"

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


process ConvertVcfToBgen {
    label "medium2"

    input:
        tuple val(chr), path(vcf_file)

    output:
        tuple path("chr*bgen"), path("chr*sample"), emit: bgen_ch
    
    script:
    """
        ${projectDir}/tools/plink --vcf $vcf_file dosage=DS --export bgen-1.2 ref-first --out chr${chr} --chr $chr
    """
}

process ConvertVcfToPlink {
    label "medium1"

    input:
        tuple val(chr), path(vcf_file)

    output:
        tuple path("chr*bed"), path("chr*bim"), path("chr*fam"), emit: bfile_per_chr_ch
    
    script:
    """    
        ${projectDir}/tools/plink2  \
        --vcf $vcf_file \
        --make-bed --out chr${chr} \
        --chr $chr  \
        --const-fid \
        --maf 0.05 --hwe 1e-06 --geno 0.05 --mac 10 \
        --extract-if-info "R2 > 0.4"

    """
}

process MergePlinkPerChr {
    label "medium2"

    //echo true
    input:
        path(plink_files)
    output:
        tuple path("merged.bed"), path("merged.bim"), path("merged.fam"), emit: bfile_ch
    shell:
    '''
    plink_exe=!{projectDir}/tools/plink
    for f in chr*bed
    do
        echo ${f%.bed} >> filelist.txt
    done

    #cat filelist.txt

    ${plink_exe} --merge-list filelist.txt --make-bed --out merged
    '''
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
	signature_matrix = "$projectDir/data/signature_matrices/" + signature_matrix_name  + "_ensg.txt.gz"
    
    if (deconvolution_method == "NA") {
        rnaquality_ch = CalculateRNAQualityScore(normalized_expression_data).view()

        CombineCovariatesRNAqual(covariates,Channel.fromPath("NA"), genotype_pcs, gte, rnaquality_ch)
        covariates_ch = CombineCovariatesRNAqual.out.covariates_ch
    } else if (deconvolution_method == "lab"){
        if (exp_type == "RNAseq" || exp_type == "RNAseq_HGNC") {
            cell_counts_ch = Channel.fromPath(params.lab_cell_perc)
            rnaquality_ch = CalculateRNAQualityScore(normalized_expression_data)
            CombineCovariatesRNAqual(covariates, cell_counts_ch, genotype_pcs, gte, rnaquality_ch)
            covariates_ch = CombineCovariatesRNAqual.out.covariates_ch
        } else {
            cell_counts_ch = Channel.fromPath(params.lab_cell_perc)
            covariates_ch = CombineCovariates(covariates,cell_counts_ch, genotype_pcs, gte)
        }
    } else {
        if (exp_type == "RNAseq" || exp_type == "RNAseq_HGNC") {
            cell_counts_ch = Deconvolution(TPM(raw_expression_data, gene_lengths), signature_matrix, deconvolution_method, exp_type)
            rnaquality_ch = CalculateRNAQualityScore(normalized_expression_data)
            CombineCovariatesRNAqual(covariates,cell_counts_ch, genotype_pcs, gte, rnaquality_ch)
            covariates_ch = CombineCovariatesRNAqual.out.covariates_ch
        } else {
            cell_counts_ch = Deconvolution(normalized_expression_data, signature_matrix, deconvolution_method, exp_type)
            covariates_ch = CombineCovariates(covariates,cell_counts_ch, genotype_pcs, gte)
        }
    }   
    
    emit:
        covariates_ch

}



