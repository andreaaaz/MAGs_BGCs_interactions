#!/usr/bin/env nextflow

params.temps = ['global','low','mid','high']
params.indir = 'metadata/'
params.outdir = 'exp/perm_temps/'
params.microbial_lineage = "mOTUs_Species_Cluster"
params.bgc_groups = "gcf"
params.nperms = 10

process PERMUTE_SITES {

  tag "$temp"

  publishDir "${params.outdir}/sites/${temp}", mode: 'copy'

  input:
  val temp
  
  output:
  tuple val(temp), path("perm_sites_${temp}_*.rds"), emit: perms
  tuple val(temp), path("real_sites_${temp}.rds"), emit: real

  script:
  """
  Rscript ${projectDir}/create_perms.r \
    --temp ${temp} \
    --indir ${params.indir} \
    --outdir ./ \
    --perms ${params.nperms} \
    --workdir ${projectDir}/ 
  """
}


process INTERACTIONS {

  tag "${temp}:${perm_file.baseName}"
  publishDir "${params.outdir}/cases/${temp}", mode: 'copy'

  input:
  tuple val(temp), path(perm_file)

  output:
  path "all_cases_*.csv"
  path "oc_filt_*.csv"

  script:
  """
  Rscript ${projectDir}/temperatures.r \
    --temp ${temp} \
    --indir ${params.indir} \
    --outdir ./ \
    --perm ${perm_file} \
    --workdir ${projectDir}/
  """
}

workflow {

    temp_ch = Channel.fromList(params.temps)

    perms = PERMUTE_SITES(temp_ch)

    perms.perms.view()
}


