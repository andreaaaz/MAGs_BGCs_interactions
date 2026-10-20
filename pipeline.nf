#!/usr/bin/env nextflow

params.indir = 'metadata/'
params.outdir = 'exp/2026-interactions/'

params.method = "binomial"
params.microbial_lineage = "mOTUs_Species_Cluster"
params.bgc_groups = "gcf"
params.quality = "8"

params.temps = ['global','low','mid','high']
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
    tuple val(temp), path("nodes_mm.csv"), emit: nodes_mm

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
    tuple val(temp), path("edges_mb.csv"),  emit: edges_mb
    tuple val(temp), path("nodes_mb.csv"),  emit: nodes_mb
    tuple val(temp), path("edges_mbm.csv"), emit: edges_mbm
    tuple val(temp), path("edges_mmr.csv"), emit: edges_mmr
    tuple val(temp), path("nodes_mmr.csv"), emit: nodes_mmr


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

    publishDir "${params.outdir}/${params.microbial_lineage}_${params.bgc_groups}/${temp}/", mode: 'copy'

    input: 
    tuple val(temp), path(mm_file), path(mmr_file)

    output: 
    tuple val(temp), path("edges_cons.csv"), emit: edges_cons
    tuple val(temp), path("nodes_cons.csv"), emit: nodes_cons

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


process TYPES_STATS {

    tag "$temp"

    publishDir "${params.outdir}/${params.microbial_lineage}_${params.bgc_groups}/type_stats/", mode: 'copy'

    input:
    tuple val(temp), 
        path(edges_mb), path(nodes_mb),
        path(edges_mmr), path(nodes_mmr),
        path(edges_mm), path(nodes_mm),
        path(edges_cons), path(nodes_cons)

    output:
    path("${temp}.png"), emit: plots

    script:
    """
    Rscript ${projectDir}/types_stats.r \
        -t ${temp} \
        --edges_mb ${edges_mb} \
        --nodes_mb ${nodes_mb} \
        --edges_mmr ${edges_mmr} \
        --nodes_mmr ${nodes_mmr} \
        --edges_mm ${edges_mm} \
        --nodes_mm ${nodes_mm} \
        --edges_cons ${edges_cons} \
        --nodes_cons ${nodes_cons} \
        -o ./${temp}.png
    """
} 

/*process TEMPS {
    tag "$types"

    publishDir "${params.outdir}/${params.microbial_lineage}_${params.bgc_groups}/temp_stats", mode: 'copy'

    input:
    
    output:
    
    script:
    """
    Rscript ${projectDir}/temps_stats.R \
        --type $\

    """
} */


workflow {
    
    temp_ch = Channel.fromList(params.temps)
    
    // interactions
    mag_bgc_out = MAG_BGC(temp_ch)
    mag_mag_out = MAG_MAG(temp_ch)
    
    // networks
    networks_mb_out = NETWORKS_MB(mag_bgc_out.oc_filt)
    networks_mm_out = NETWORKS_MM(mag_mag_out.oc_filt_mm)

    // consensus
    consensus_inputs = networks_mm_out.edges_mm.join(
        networks_mb_out.edges_mmr
    )
    consensus_out = CONSENSUS(consensus_inputs)

    //type stats
    // group edges and nodes of each network type 
    mb_files = networks_mb_out.edges_mb.join(networks_mb_out.nodes_mb)
    mmr_files = networks_mb_out.edges_mmr.join(networks_mb_out.nodes_mmr)
    mm_files = networks_mm_out.edges_mm.join(networks_mm_out.nodes_mm)
    cons_files = consensus_out.edges_cons.join(consensus_out.nodes_cons)

    // join them
    type_stats_inputs = mb_files
        .join(mmr_files)
        .join(mm_files)
        .join(cons_files)

    TYPES_STATS(type_stats_inputs)

}


