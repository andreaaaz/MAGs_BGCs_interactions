#!/usr/bin/env nextflow

params.indir = 'metadata/'
params.outdir = 'exp/2026-interactions/'

params.method = "binomial"
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
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/interactions_MAG-BGC.r \
        -m ${params.microbial_lineage} \
        -b ${params.bgc_groups} \
        -t ${temp} \
        -i ${params.indir} \
        -o ./ \
        -w ${projectDir}/ \
        -e ${params.method}
    """
}


process MAG_MAG {
        
    tag "$temp"
        
    publishDir "${params.outdir}/${params.microbial_lineage}/${temp}/", mode: 'copy'
    
    input:
    val temp

    output:
    tuple val(temp), path("oc_filt.csv"), emit: oc_filt_mm
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/interactions_MAG-MAG.r \
        -m ${params.microbial_lineage} \
        -t ${temp} \
        -i ${params.indir} \
        -o ./ \
        -w ${projectDir}/ \
        -e ${params.method}
    """
}   

process NETWORKS_MM {

    tag "$temp"

    publishDir "${params.outdir}/${params.microbial_lineage}/${temp}/", mode: 'copy'

    input:
    tuple val(temp), path(oc_file)

    output:
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/networks_mm.r \
        -m ${params.microbial_lineage} \
        -i ${params.indir} \
        -f ${oc_file} \
        -o ./ 
    """

}

process NETWORKS_MB {
    
    tag "$temp"

    publishDir "${params.outdir}/${params.microbial_lineage}_${params.bgc_groups}/${temp}/", mode: 'copy'

    input:
    tuple val(temp), path(oc_mb)

    output:
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/networks_mb.r \
        -m ${params.microbial_lineage} \
        -b ${params.bgc_groups} \
        -i ${params.indir} \
        -o ./ \
        -w ${projectDir}/ \
        -f ${oc_mb} 
    """

}


workflow {
    
    temp_ch = Channel.fromList(params.temps)
    
    // interacciones 
    mag_bgc_out = MAG_BGC(temp_ch)
    mag_mag_out = MAG_MAG(temp_ch)
    
    // redes
    NETWORKS_MB(mag_bgc_out.oc_filt)
    NETWORKS_MM(mag_mag_out.oc_filt_mm)

}

