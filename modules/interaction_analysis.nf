#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * run QTL mapping
 */
process ieQTL_mapping {
    tag "Chunk: $chunk"

    input:
    tuple path(tmm_expression), path(bed), path(bim), path(fam), path(covariates), path(limix_annotation), path(gte), path(qtls_to_test), val(covariate_to_test), val(chunk)
    

    output:
    //path "limix_out/*"

    shell:
    '''
     geno=!{bed}
     plink_base=${geno%.bed}
     outdir=${PWD}/limix_out/
     mkdir $outdir
     ls $PWD

     python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --plink ${plink_base} \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{tmm_expression} \
      -smf !{gte} \
      -fvf !{qtls_to_test} \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -gr !{chunk} \
      -np 0 \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001
      
      echo !{params.outdir}
      if [ ! -d !{params.outdir}/limix_output/ ]
      then
      	mkdir !{params.outdir}/limix_output/
      fi
      if [ ! -z "$(ls -A ${outdir}/)" ]
      then
      	cp ${outdir}/* !{params.outdir}/limix_output/
      else
	echo "No limix output to copy"      
      fi
      ls -la ${outdir}/  
      
    '''
}


