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
    --norm_expfile exp_data_preprocessed.txt.gz \
    --gte LLD_gte.txt \
    --exp_platform RNAseq \
    --cohort_name LLD \
    --covariates LLD_covariates.txt \
    --covariate_to_test gender_F1M2 \
    --genes_to_test signif_eqtl_genes.txt \
    --outdir LLD_interaction_res \
    --plink2_executable plink2 \
    -profile slurm \
    -resume

  Mandatory arguments:
    --cohort_name                 Name of the cohort.
    --bfile                       Path to imputed genotype files as used in eQTLGen main analysis in plink bed/bim/fam format (without extensions bed/bim/fam).
    --raw_expfile                 Path to the un-preprocessed gene expression matrix (genes/probes in the rows, samples in the columns). Can be from RNA-seq experiment or from array. NB! For Affymetrix arrays (AffyU219, AffyExon) we assume that standard preprocessing and normalisation is already done.
    --norm_expfile                Path to the normalized gene expression matrix produced by DataQC step in the primary eQTLGen pipeline.
    --gte                         Genotype-to-expression linking file. Tab-delimited, no header. First column: sample ID for genotype data. Second column: corresponding sample ID for gene expression data. Can be used to filter samples from the analysis.
    --exp_platform                Expression platform. Should be one of HT12v3, HT12v4, HuRef8, RNAseq, AffyU219, AffyHumanExon, RNAseq_HGNC. 
    --covariates                  Path to the covariate file containing the covariate to test for interaction, sex, age and cohort-specific covariates. Tab-delimited, with header (Samples in rows, covariates in columns). NB! Covariate ids should be the same as used in the genotype data.
    --covariate_to_test           Name of covariate to test for interaction. Should be one of the column names of the covariate file.
    --outdir                      Path to the output directory.

  Optional arguments
    --signature_matrix_name       Name of the signature matrix using in cell type deconvolution. Either LM22 or ABIS
    --deconvolution_method        Name of the cell type deconvolution method. Either nnls or dtangle
    --plink2                      Plink2 executable
    
  """.stripIndent()
}

/*
 * Parameters
 */

params.bfile = ''
params.vcf_dir = ''
params.bgen_dir = ''
//params.bfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output//chr2"
//params.vcf_dir = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/input/postimpute/"
//params.bgen_dir = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/"

params.genes_to_test = ''
params.qtls_to_test = ''
params.preadjust = false
params.num_expr_PCs = 25
params.signature_matrix_name = "LM22"
params.deconvolution_method = "dtangle"
params.num_perm = 0

params.run_stratified = false
params.plink2_executable = "plink2"
/*
 * Channel declarations
 */

raw_expr_ch = Channel.fromPath(params.raw_expfile)
filt_exp_ch = Channel.fromPath(params.norm_expfile)
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)
gte_ch = Channel.fromPath(params.gte)

annotation_ch = Channel.fromPath("$projectDir/data/LimixAnnotationFile.GRCh38.110.txt.gz")
gene_lengths_ch = Channel.fromPath("$projectDir/data/GeneLengths_GRCh38.110_ensg.txt.gz")

//Channel.fromPath(params.genes_to_test).set { genes_to_test_ch } 

Channel
    .fromPath(params.chunk_file)
    .splitCsv( header: false )
    .map { row -> tuple(row[0].split(':')[0], row[0]) }
    .set { chunk_ch }


if (params.bgen_dir != '') {
  Channel.from(1..22)
  .map { chr -> tuple("$chr", file("${params.bgen_dir}/*chr${chr}.bgen"), file("${params.bgen_dir}/*chr${chr}.sample")) }
  .ifEmpty { exit 1, "Input .bgen files not found!" }
  .set { chr_bgen_pairs }
} else if (params.vcf_dir != '') {
  Channel.from(1..22)
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


include { PREPARE_COVARIATES; PrepareAnnotation; NormalizeExpression; ConvertVcfToBgen; ConvertVcfToPlink; MergePlinkPerChr } from './modules/prepare_data.nf'
//include { IeQTLmappingPerGeneTMP; IeQTLmappingPerSNPGene; IeQTLmappingPerGene; IeQTLmappingPerGeneNoChunks; IeQTLmappingPerGeneBgen; FilterGenesToTest } from './modules/interaction_analysis.nf'
include { RUN_INTERACTION_QTL_MAPPING; IeQTLmapping; IeQTLmapping_InteractionCovariates; SplitCovariates; PreadjustExpression } from './modules/interaction_analysis2.nf'
include { RUN_STRATIFIED_ANALYSIS; RunEqtlMappingPerGenePlink } from './modules/stratified_analysis.nf'

/* 
 * Analysis
 */

workflow {
    /*
     * Prepare expression and covariate data
     */
    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
    
    NormalizeExpression(raw_expr_ch, filt_exp_ch, params.exp_platform, gte_ch )
    norm_exp_ch = NormalizeExpression.out.norm_expression_table
    covariates_ch = PREPARE_COVARIATES(params.exp_platform, raw_expr_ch, norm_exp_ch, params.signature_matrix_name, params.deconvolution_method,covars_ch, gene_lengths_ch, annotation_ch, params.genotype_pcs, params.gte)
    //norm_exp_ch = filt_exp_ch
    //covariates_ch = Channel.fromPath("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/covariates.combined.txt")
    
    /*
     * Prepare genotype data
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

      Channel.empty()
        .set { bfile_ch }

      if (params.bfile == '') {
        ConvertVcfToPlink(chr_vcf_pairs)
        MergePlinkPerChr(ConvertVcfToPlink.out.bfile_per_chr_ch.collect()).set {bfile_ch}

        //interaction_ch = norm_exp_ch.combine(MergePlinkPerChr.out.bfile_ch).combine(covariates_ch).combine(annotation_ch).combine(gte_ch).combine(Channel.fromPath(params.genes_to_test)).combine(Channel.of(params.covariate_to_test)).combine(chunk_ch.map { it[1] })
      } else {
        Channel
          .from(params.bfile)
          .ifEmpty { exit 1, "Input plink prefix not found!" }
          .map { genotypes -> [file("${genotypes}.bed"), file("${genotypes}.bim"), file("${genotypes}.fam")]}
          .set { bfile_ch }
      }

      //bfile_ch.view()
      /*
       * Run interaction analysis
       */      
      
      //interaction_ch = norm_exp_ch.combine(bfile_ch).combine(covariates_ch).combine(annotation_ch).combine(gte_ch).combine(genes_to_test_ch).combine(Channel.of(params.covariate_to_test)).combine(chunk_ch.map { it[1] })
      //results_ch = IeQTLmappingPerGene(interaction_ch)
      RUN_INTERACTION_QTL_MAPPING(norm_exp_ch, bfile_ch, covariates_ch, annotation_ch, Channel.of(params.covariate_to_test), chunk_ch.map { it[1] })

      /*
       * Stratified analysis 
       */
       if (params.run_stratified){
        RUN_STRATIFIED_ANALYSIS(norm_exp_ch, bfile_ch, covariates_ch, annotation_ch, gte_ch, genes_to_test_ch, chunk_ch)
       }
    //}
    
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
