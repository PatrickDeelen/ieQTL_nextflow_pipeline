#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

limix_path="singularity exec --bind /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/ /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_limix/limixAug23.simg python /limix_qtl/Limix_QTL/"

/*
 * run QTL mapping
 */
process ieQTL_mapping {

    input:
    path tmm_expression
    path covariates
    path annotation
    path interactions_outdir
    val covariate_name
    val num_perm
    

    output:
    

    script:
    """
     ls $interactions_outdir

     ${limix_path}/run_interaction_QTL_analysis.py \
     --plink params.bfile \
      -af $annotation \
      -cf $covariates \
      -pf $tmm_expression \
      -smf params.gte \
      -fvf params.qtls_to_test \
      -od $interactions_outdir \
      --interaction_term ${covariate_name} \
      -np $num_perm \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001
      
      ${limix_path}/post_processing/minimal_interaction_postprocess.py \
      -id $interactions_outdir \
      -od $interactions_outdir \
      -sfo
      
      
    """
}

