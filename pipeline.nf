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
    path "*filt.csv"

    script:
    """
    Rscript ${projectDir}/interactions_MAG-BGC.r \
        -m ${params.microbial_lineage} \
        -b ${params.bgc_groups} \
        -t ${temp} \
        -i ${params.indir} \
        -o .
    """
}


process MAG_MAG {
        
    tag "$temp"
        
    publishDir "${params.outdir}/${params.microbial_lineage}_${params.microbial_lineage}/${temp}/", mode: 'copy'
    
    input:
    val temp

    output:
    path "*filt.csv"

    script:
    """
    Rscript ${projectDir}/interactions_MAG-MAG.r \
        -m ${params.microbial_lineage} \
        -t ${temp} \
        -i ${params.indir} \
        -o .
    """
}   


workflow {
    
    temp_ch = Channel.from(params.temps)
    
    MAG_BGC(temp_ch)
    MAG_MAG(temp_ch)

}
