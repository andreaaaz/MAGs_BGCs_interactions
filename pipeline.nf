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
    tuple val(temp), path("edges_mm.csv"), emit: edges_mm
    path("nodes_mm.csv"), emit: nodes_mm

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
    path("nodes_mb.csv"), emit: nodes_mb
    path("edges_mb.csv"), emit: edges_mb
    path("edges_mbm.csv"), emit: edges_mbm
    tuple val(temp), path("edges_mmr.csv"), emit: edges_mmr
    path("nodes_mmr.csv"), emit: nodes_mmr

    script:
    """
    Rscript ${projectDir}/networks_mb.r \
        -m ${params.microbial_lineage} \
        -b ${params.bgc_groups} \
        -t ${temp} \
        -q ${params.quality} \
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
    tuple val(temp), path(mm_file), path(mmr_file)

    output: 
    tuple val(temp), path("edges_cons.csv"), emit: edges_cons
    path("nodes_cons.csv"), emit: nodes_cons

    script:
    """
    Rscript ${projectDir}/consensus.R \
        -t ${temp} \
        -i ${params.indir} \
        -m ${mm_file} \
        -r ${mmr_file} \
        -o ./ \
        -w ${projectDir}/
    """
}

/* process TEMPS_STATS {

    tag "$temp"

    publishDir "${params.outdir}/stats/", mode: 'copy'

    input:
    val temp // necesito los nodes y edges de las mm, las mb y las mmr y las consensus de cada temperatura 

    output:
    // deben ser figuras

    script:
    """
    Rscript ${projectDir}/temps_stats.r \
        -t ${temp} \
        -i ${params.indir} \
        -o ./   
// sigue en construccion
    """
} */

/* process TYPES_STATS {

    tag "$temp"

    publishDir "${params.outdir}/stats/", mode: 'copy'

    input:
    val temp // necesito los nodes y edges de las mm, las mb y las mmr y las consensus de cada temperatura 

    output:
    // deben ser figuras

    script:
    """
    Rscript ${projectDir}/types_stats.r \
        -t ${temp} \
        -i ${params.indir} \
        -o ./ 
    // sigue en contruccion
    """
} */


workflow {
    
    temp_ch = Channel.fromList(params.temps)
    
    // interacciones 
    mag_bgc_out = MAG_BGC(temp_ch)
    mag_mag_out = MAG_MAG(temp_ch)
    
    // redes
    networks_mb_out = NETWORKS_MB(mag_bgc_out.oc_filt)
    networks_mm_out = NETWORKS_MM(mag_mag_out.oc_filt_mm)

    // consensus
    consensus_inputs = networks_mm_out.edges_mm.join(
        networks_mb_out.edges_mmr
    )
    consensus_out = CONSENSUS(consensus_inputs)

    //stats
}


