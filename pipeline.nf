#!/usr/bin/env nextflow

params.indir = 'metadata/'
params.outdir = 'exp/2026-interactions/'

params.microbial_lineage = "mOTUs_Species_Cluster"
params.bgc_groups = "gcf"

params.temps = ['global','low','mid','high']


process MAG_BGC {

    tag "$temp"

    publishDir "${params.outdir}/${params.microbial_lineage}_${params.bgc_groups}/${temp}/", mode: 'copy'

    input:
    val temp

    output:
    tuple val(temp), path("oc_filt.csv"), emit: oc_filt
    path "*.csv", emit: all_csvs

    script:
    """
    Rscript ${projectDir}/MAGs_BGCs_interactions/interactions_MAG-BGC.r \
        -m ${params.microbial_lineage} \
        -b ${params.bgc_groups} \
        -t ${temp} \
        -i ${params.indir} \
        -o ./
    """
}


process MAG_MAG {
        
    tag "$temp"
        
    publishDir "${params.outdir}/${params.microbial_lineage}_${params.microbial_lineage}/${temp}/", mode: 'copy'
    
    input:
    val temp

    output:
    tuple 
 

    script:
    """
    Rscript ${projectDir}/MAGs_BGCs_interactions/interactions_MAG-MAG.r \
        -m ${params.microbial_lineage} \
        -t ${temp} \
        -i ${params.indir} \
        -o ./
    """
}   

process NETWORKS_MB {
    
    tag "$temp"

    publishDir "${params.outdir}/${params.microbial_lineage}_${params.bgc_groups}/${temp}/", mode: 'copy'

    input:
    tuple val(temp), path(oc_file)

    output:
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/MAGs_BGCs_interactions/networks_mb.r \
        -m ${params.microbial_lineage} \
        -b ${params.bgc_groups} \
        -f ${oc_file} \
        -o ./
    """

}


workflow {
    
    temp_ch = Channel.of(params.temps)
    
    mag_bgc_out = MAG_BGC(temp_ch)
    
    NETWORKS_MB(mag_bgc_out.oc_filt)

    MAG_MAG(temp_ch)

}

