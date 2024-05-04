#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * run QTL mapping per SNP-Gene pair or for all SNPs around the gene depending on command line parameters 
 */
process IeQTLmapping {
    tag "Chunk: $chunk"
    //echo true

    input:
    tuple path(tmm_expression), path(bed), path(bim), path(fam), path(covariates), path(limix_annotation), val(covariate_to_test), val(chunk), path(qtl_ch)
    

    output:
    path "limix_out/*"

    shell:
    '''
    geno=!{bed}
    plink_base=${geno%.bed}
    outdir=${PWD}/limix_out/
    mkdir -p $outdir

    awk 'BEGIN {OFS="\\t"}; {print $2, $2}' !{fam} > gte.txt

    qtls=!{params.qtls_to_test}
    genes=!{params.genes_to_test}
    if [ "${#qtls}" -gt 1 ]
    then 
        arg_line="-fvf !{qtl_ch}"
    else
        arg_line="-ff !{qtl_ch}"
    fi

    #python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \

    #eval "$(conda shell.bash hook)"
    #conda activate py39
    #python /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_limix/limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
    
    python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --plink ${plink_base} \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{tmm_expression} \
      -smf gte.txt \
      ${arg_line} \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -gr !{chunk} \
      -np !{params.num_perm} \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001 \
      --write_permutations --write_zscore
      
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
      
    '''
}


/*
 * run QTL mapping per SNP-Gene pair or for all SNPs around the gene depending on command line parameters 
 */
process IeQTLmapping_InteractionCovariates {
    tag "Chunk: $chunk"

    input:
    tuple path(preadjusted_expression), path(bed), path(bim), path(fam), path(interaction_covariates), path(limix_annotation), val(covariate_to_test), val(chunk)
    

    output:
    //path "limix_out/*"

    shell:
    '''
    geno=!{bed}
    plink_base=${geno%.bed}
    outdir=${PWD}/limix_out/
    mkdir $outdir

    awk 'BEGIN {OFS="\\t"}; {print $2, $2}' !{fam} > gte.txt

    qtls=!{params.qtls_to_test}
    genes=!{params.genes_to_test}
    if [ "${#qtls}" -gt 1 ]
    then 
        arg_line="-fvf $qtls_ch"
    else
        arg_line="-ff $qtls_ch"
    fi

    #python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \

    #eval "$(conda shell.bash hook)"
    #conda activate py39
    #python /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_limix/limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
    python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --plink ${plink_base} \
      -af !{limix_annotation} \
      -cf !{interaction_covariates} \
      -pf !{preadjusted_expression} \
      -smf gte.txt \
      $arg_line \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -gr !{chunk} \
      -np !{params.num_perm} \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001 \
      --interaction_covariates \
      --write_permutations --write_zscore
      
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


/*
 * split covariates into those that should be regressed out before the interaction analysis and those that will be included in the interaction model together with their interaction with genotypes
 */
process SplitCovariates {
    tag "Split covariates"
    label "short"
    //publishDir params.outdir, mode: 'copy'

    input:
    path covariates_path

    output:
    path ("covariates_preadjust.txt"), emit: linear_covariates_ch
    path ("covariates_interaction.txt"), emit: interaction_covariates_ch

    shell:
    '''
    if [ !{params.exp_platform} = "RNAseq" -o !{params.exp_platform} = "RNAseq_HGNC" ]; then
          interaction_covariates=`echo "AvgExprCorrelation,!{params.covariate_to_test}"`
    else 
          interaction_covariates = `echo "!{params.covariate_to_test}"`
    fi
    if [ !{params.cell_perc_interactions} = false ]; then
        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates} covariates_interaction.txt
        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates} covariates_preadjust.txt -v
    else 
        cell_types=`head -1  ../../run1/results/RS/cell_counts.txt | cut -f2- | sed "s:\t:,:g"`
        interaction_covariates2=`echo "$interaction_covariates,$cell_types"`

        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates2} covariates_interaction.txt
        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates2} covariates_preadjust.txt -v

    fi
    '''
}

/*
 * split covariates into those that should be regressed out before the interaction analysis and those that will be included in the interaction model together with their interaction with genotypes
 */
process PreadjustExpression {
    tag "Preadjust expression"
    label "short"
    //publishDir params.outdir, mode: 'copy'

    input:
    path expression_path
    path covariates_path
    path pcs_path

    output:
    path ("expression_corrected.noINT.txt")

    shell:
    '''
    n=!{params.num_expr_PCs}
    n2=$((n + 1))
    cut -f1-$n2 !{pcs_path} > pcs.txt
    Rscript !{x}/bin/regress_linear_covariates.R !{expression_path} !{covariates_path} ./ pcs.txt
    '''
}

process ConvertIeQTLsToText {
    echo true

    input:
    path limix_out_files
    
    output:
    

    script:
    """
    mkdir limix_out_text/
    cp $limix_out_files limix_out_text/
    python /limix_qtl/Limix_QTL/post_processing/minimal_interaction_postprocess.py \
      -id limix_out_text/ \
      -od  limix_out_text/ \
      -sfo
    
    gzip limix_out_text/*txt 
    mv limix_out_text/*txt.gz ${params.outdir}/limix_output/

    """
}

process GetTopIeQTLs {
    echo true

    input:
    path limix_folder
    
    output:
    

    script:
    """
    """
}
workflow RUN_INTERACTION_QTL_MAPPING {   
    take:
        tmm_expression
	    plink_geno
        covariates_ch
	    limix_annotation
        covariate_to_test
	    chunk
        qtl_ch

    main:

    
    if (params.preadjust){            
        SplitCovariates(covariates_ch)
        PreadjustExpression(tmm_expression, SplitCovariates.out.linear_covariates_ch, params.expr_pcs)

        interaction_ch = PreadjustExpression.out.combine(plink_geno).combine(SplitCovariates.out.interaction_covariates_ch).combine(limix_annotation).combine(covariate_to_test).combine(chunk)
        ConvertIeQTLsToText(IeQTLmapping_InteractionCovariates(interaction_ch).collect())
    } else {
        interaction_ch = tmm_expression.combine(plink_geno).combine(covariates_ch).combine(limix_annotation).combine(covariate_to_test).combine(chunk).combine(qtl_ch)
        ConvertIeQTLsToText(IeQTLmapping(interaction_ch).collect())
    
    }   
    
    //emit:
    //    results_ch

}