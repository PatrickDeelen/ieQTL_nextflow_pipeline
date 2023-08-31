#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


/*
 * run QTL mapping
 */
process ieQTL_mapping {
    input:
    path tmm_expression
    path genotypes
    path covariates
    path annotation
    path interactions_outdir
    path gte
    path qtls_to_test
    val covariate_name
    val num_perm
    

    output:
    

    script:
    """
     ls \$PWD
     echo \$PWD/$covariates
     python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --plink $genotypes \
      -af \$PWD/$annotation \
      -cf \$PWD/$covariates \
      -pf \$PWD/$tmm_expression \
      -smf \$PWD/$gte \
      -fvf \$PWD/$qtls_to_test \
      -od \$PWD/$interactions_outdir \
      --interaction_term ${covariate_name} \
      -np $num_perm \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001
      
      
      
    """
}

