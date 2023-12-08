#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
    log.info"""
    =======================================================
     Interaction analysis v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run InteractionAnalysis.nf \
        --bfile LLD_genotypes_flt \
        --raw_expfile LLD_expression_data.txt.gz \
        --norm_expfile LLD_expression_data.TMM.txt.gz \
        --gte LLD_gte.txt \
        --exp_type rnaseq \
        --cohort_name LLD \
        --covariates LLD_covariates.txt \
        --covariate_to_test gender_F1M2 \
        --gtf Homo_sapiens.GRCh37.75.gtf.gz \
        --genes_to_test signif_eqtl_genes.txt \
        --genome_build GRCh37 \
        --outdir LLD_interaction_res \
	--plink2_executable plink2 \
        -profile slurm \
        -resume

    Mandatory arguments:
      --cohort_name                 Name of the cohort.
      --genome_build                Genome build of the cohort. Either hg18, GRCh36, hg19, GRCh37, hg38 or GRCh38.
      --bfile                       Path to imputed genotype files as used in eQTLGen main analysis in plink bed/bim/fam format (without extensions bed/bim/fam).
      --raw_expfile                 Path to the un-preprocessed gene expression matrix (genes/probes in the rows, samples in the columns). Can be from RNA-seq experiment or from array. NB! For Affymetrix arrays (AffyU219, AffyExon) we assume that standard preprocessing and normalisation is already done.
      --norm_expfile                Path to the normalized gene expression matrix before correcting for PCs. E.g. TMM-transformed matrix for RNA-seq data.
      --gte                         Genotype-to-expression linking file. Tab-delimited, no header. First column: sample ID for genotype data. Second column: corresponding sample ID for gene expression data. Can be used to filter samples from the analysis.
      --exp_type                    Expression estimation method. Either rnaseq or array.
      --covariates                  Path to the covariate file containing the covariate to test for interaction, sex, age and cohort-specific covariates. Tab-delimited, with header (Samples in rows, covariates in columns). NB! Covariate ids should be the same as used in the expression table.
      --covariate_to_test           Name of covariate to test for interaction. Should be one of the column names of the covariate file.
      --outdir                      Path to the output directory.
      --gtf                         Path to the gtf annotation file of the same build as used in the expression matrix
      --outdir                      Path to the output folder

    Optional arguments
      --signature_matrix_name       Name of the signature matrix using in cell type deconvolution. Either LM22 or ABIS
      --deconvolution_method        Name of the cell type deconvolution method. Either nnls or dtangle
      
    """.stripIndent()
}

//params.raw_expfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/expression_data.all.txt.gz"
params.raw_expfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/expression_data.txt"
params.norm_expfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/LLD.gene_read_counts_BIOS_and_LLD_passQC.TMM.txt.gz"
params.bfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input//LLD_genotypes_flt"
params.covariates = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/v2/LLD_age_sex.txt"

params.gte = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/LLD_gte.txt"
params.outdir = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/"
params.covariate_to_test = "gender_F1M2"
params.gtf = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/Homo_sapiens.GRCh37.75.gtf.gz"
params.genes_to_test = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/genes_to_test.txt"
params.genotype_pcs = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/LLD_PCs.txt"

params.genome_build = "GRCh37"
params.exp_type = "rnaseq"
params.cohort_name = "LLD"
params.signature_matrix_name = "LM22"
params.deconvolution_method = "nnls"
params.num_perm = 0

params.plink2_executable = "plink2"

raw_expr_ch = Channel.fromPath(params.raw_expfile)
norm_exp_ch = Channel.fromPath(params.norm_expfile)
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)
gtf_annotation_ch = Channel.fromPath(params.gtf)

annotation_ch = Channel.fromPath("$projectDir/data/LimixAnnotationFile.txt.gz")
gene_lengths_ch = Channel.fromPath("$projectDir/data/GeneLengths.txt.gz")

Channel
    .from("$projectDir/data/ChunkingFile.txt")
    .splitCsv( header: false )
    .map { row -> row[0] }
    .set { chunk_ch }

Channel
    .from(params.bfile)
    .ifEmpty { exit 1, "Input plink prefix not found!" }
    .map { genotypes -> [file("${genotypes}.bed"), file("${genotypes}.bim"), file("${genotypes}.fam")]}
    .set { bfile_ch }


include { TMM_TRANSFORM_EXPRESSION; PREPARE_COVARIATES; PrepareAnnotation } from './modules/prepare_data.nf'
include { IeQTLmappingPerSNPGene; IeQTLmappingPerGene } from './modules/interaction_analysis.nf'

workflow {
    //norm_exp_ch = TMM_TRANSFORM_EXPRESSION(raw_expr_ch)

    //PrepareAnnotation(gtf_annotation_ch)

    covariates_ch = PREPARE_COVARIATES(params.exp_type, raw_expr_ch, norm_exp_ch, params.signature_matrix_name, params.deconvolution_method,covars_ch, gene_lengths_ch, annotation_ch, params.genotype_pcs, params.gte)
    
    //to run in chunks:
    //PrepareAnnotation.out.chunks_ch.splitCsv( header: false ).map { row -> row[0] }.set { chunk_ch }
    eqtl_ch = norm_exp_ch.combine(bfile_ch).combine(covariates_ch).combine(annotation_ch).combine(Channel.fromPath(params.gte)).combine(Channel.fromPath(params.genes_to_test)).combine(Channel.of(params.covariate_to_test)).combine(chunk_ch)
    results_ch = IeQTLmappingPerGene(eqtl_ch)

    //without chunks:
    //eqtl_ch = norm_exp_ch.combine(bfile_ch).combine(covariates_ch).combine(PrepareAnnotation.out.annotation_ch).combine(Channel.fromPath(params.gte)).combine(Channel.fromPath(params.qtls_to_test)).combine(Channel.of(params.covariate_to_test))
    //results_ch = ieQTL_mapping_no_chunks(eqtl_ch)

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
