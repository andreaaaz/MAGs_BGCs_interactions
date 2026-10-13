#!/usr/bin/env nextflow

params.indir = 'metadata/'
params.outdir = 'exp/2026-interactions/'

params.method = "binomial"
params.microbial_lineage = "mOTUs_Species_Cluster"
params.bgc_groups = "gcf"
params.quality = "8"

params.temps = ['global', 'low_depth','low','mid','high']
params.types = ['MAG-BGC', 'MAG-MAG', 'MAG-MAG-rec', 'consensus']

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
        -q ${params.quality} \
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
        -q ${params.quality} \
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
    tuple val(temp), path("edges_mm.csv") // necesito que los edges_mm.csv vayan a consensus 
                                        // y tambien se tiene que ir a temp_stats y types_stats
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
    tuple val(temp), path("edges_mm.csv") // necesito que los edges_mm.csv (este tambien crea mm) vayan a consensus
                                        // y tambien se tiene que ir a temp_stats y types_stats    
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/networks_mb.r \
        -m ${params.microbial_lineage} \
        -b ${params.bgc_groups} \
        -t ${temp} \
        -i ${params.indir} \
        -o ./ \
        -w ${projectDir}/ \
        -f ${oc_mb} 
    """

}

process CONSENSUS {
    tag "$temp" 

    publishDir "${params.outdir}/${params.microbial_lineage}_${params.bgc_groups}/${temp}/consensus/", mode: 'copy'

    input: 
    touple val (temp) // aqui debe deben ir los dos edges_mm.csv de MAG-BGC y MAG-MAG

    output: 
    path "*.csv" // aqui ademas tiene que ir

    script:
    """
    Rscript ${projectDir}/consensus.R \ 
        -t ${temp} \
        -i ${params.indir} \ 
        -o ./ \
        -w ${projectDir}/ 
    """
}

process TEMPS_STATS {

    tag "$temp"

    publishDir "${params.outdir}/stats/", mode: 'copy'

    input:
    val temp

    output:
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/temps_stats.r \
        -t ${temp} \
        -i ${params.indir} \
        -o ./ 
    """
}

process TYPES_STATS {

    tag "$temp"

    publishDir "${params.outdir}/stats/", mode: 'copy'

    input:
    val temp

    output:
    path "*.csv"

    script:
    """
    Rscript ${projectDir}/types_stats.r \
        -t ${temp} \
        -i ${params.indir} \
        -o ./ 

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


