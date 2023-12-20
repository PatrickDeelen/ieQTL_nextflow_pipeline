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


params.bfile = ''
//params.vcf_dir = ''
params.bgen_dir = ''
//params.bfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output//chr2_v2"
params.vcf_dir = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/postimpute/"
//params.bgen_dir = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/"


params.raw_expfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/expression_data.all.txt.gz"
//params.raw_expfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/expression_data.txt"
params.norm_expfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/exp_data_preprocessed2.txt"
params.covariates = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/v2/LLD_age_sex.txt"

params.gte = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/LLD_gte.txt"
params.outdir = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/"
params.covariate_to_test = "gender_F1M2"
params.gtf = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/Homo_sapiens.GRCh37.75.gtf.gz"
params.genes_to_test = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/genes_to_test.txt"

params.genotype_pcs = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/LLD_PCs.txt"
params.check_sex = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/SexCheck.txt"
params.geno_fam = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input//LLD_genotypes_flt.fam"

params.genome_build = "GRCh37"
params.platform = "RNAseq"
params.cohort_name = "LLD"
params.signature_matrix_name = "LM22"
params.deconvolution_method = "nnls"
params.num_perm = 0

params.plink2_executable = "plink2"

/*
 * Channel declarations
 */

raw_expr_ch = Channel.fromPath(params.raw_expfile)
filt_exp_ch = Channel.fromPath(params.norm_expfile)
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)
gtf_annotation_ch = Channel.fromPath(params.gtf)
gte_ch = Channel.fromPath(params.gte)
annotation_ch = Channel.fromPath("$projectDir/data/LimixAnnotationFile.txt.gz")
gene_lengths_ch = Channel.fromPath("$projectDir/data/GeneLengths.txt.gz")

Channel
    .fromPath("$projectDir/data/ChunkingFile_test2.txt")
    .splitCsv( header: false )
    .map { row -> tuple(row[0].split(':')[0], row[0]) }
    .set { chunk_ch }


if (params.bgen_dir != '') {
  Channel.from(2..3)
  .map { chr -> tuple("$chr", file("${params.bgen_dir}/*chr${chr}.bgen"), file("${params.bgen_dir}/*chr${chr}.sample")) }
  .ifEmpty { exit 1, "Input .bgen files not found!" }
  .set { chr_bgen_pairs }
} else if (params.vcf_dir != '') {
  Channel.from(2..3)
  .map { chr -> tuple("$chr", file("${params.vcf_dir}/*chr${chr}.filtered.vcf.gz")) }
  .ifEmpty { exit 1, "Input .vcf.gz files not found!" }
  .set { chr_vcf_pairs }
  Channel.empty()
    .set { chr_bgen_pairs }
} else {
  Channel
    .from(params.bfile)
    .ifEmpty { exit 1, "Input plink prefix not found!" }
    .map { genotypes -> [file("${genotypes}.bed"), file("${genotypes}.bim"), file("${genotypes}.fam")]}
    .set { bfile_ch }
} 


include { PREPARE_COVARIATES; PrepareAnnotation; NormalizeExpression; NormalizeExpressionV2; ConvertVcfToBgen; ConvertVcfToPlink; MergePlinkPerChr } from './modules/prepare_data.nf'
include { IeQTLmappingPerSNPGene; IeQTLmappingPerGene; IeQTLmappingPerGeneNoChunks; IeQTLmappingPerGeneBgen } from './modules/interaction_analysis.nf'

/* 
 * Analysis
 */

workflow {
    /*
     * Prepare expression and covariate data
     */
    //NormalizeExpression(raw_expr_ch, params.platform, gte_ch, Channel.fromPath(params.check_sex), Channel.fromPath(params.geno_fam) )
    NormalizeExpressionV2(raw_expr_ch, filt_exp_ch, params.platform, gte_ch )
    norm_exp_ch = NormalizeExpressionV2.out.norm_expression_table

    //covariates_ch = PREPARE_COVARIATES(params.platform, raw_expr_ch, norm_exp_ch, params.signature_matrix_name, params.deconvolution_method,covars_ch, gene_lengths_ch, annotation_ch, params.genotype_pcs, params.gte)
    covariates_ch = Channel.fromPath("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/covariates.combined.txt")
    
    /*
     * Prepare genotype data and run ieQTL mapping
     */
    
    /*if (params.bfile == '') {
      println "No plink genotypes, will use bgen"
      if (params.bgen_dir == ''){ // Genotype convertion
        println "Converting VCF to bgen"
        ConvertVcfToBgen(chr_vcf_pairs)
        chr_bgen_pairs = ConvertVcfToBgen.out.bgen_ch.view()
      }

      // ieQTL mapping
      chunk_geno_ch = chunk_ch.join(chr_bgen_pairs).view()
      eqtl_ch = norm_exp_ch.combine(covariates_ch).combine(annotation_ch).combine(gte_ch).combine(Channel.fromPath(params.genes_to_test)).combine(Channel.of(params.covariate_to_test)).combine(chunk_geno_ch)
      results_ch = IeQTLmappingPerGeneBgen(eqtl_ch)

    } else { // Plink genotypes
      */

      if (params.bfile == '') {
        ConvertVcfToPlink(chr_vcf_pairs)
        MergePlinkPerChr(ConvertVcfToPlink.out.bfile_per_chr_ch.collect().view())
      }
      eqtl_ch = norm_exp_ch.combine(MergePlinkPerChr.out.bfile_ch).combine(covariates_ch).combine(annotation_ch).combine(gte_ch).combine(Channel.fromPath(params.genes_to_test)).combine(Channel.of(params.covariate_to_test)).combine(chunk_ch.map { it[1] })
      results_ch = IeQTLmappingPerGene(eqtl_ch)

      //run without chunks:
      //eqtl_ch = norm_exp_ch.combine(bfile_ch).combine(covariates_ch).combine(annotation_ch).combine(gte_ch).combine(Channel.fromPath(params.genes_to_test)).combine(Channel.of(params.covariate_to_test))
      //results_ch = IeQTLmappingPerGeneNoChunks(eqtl_ch)

    //}
    
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
