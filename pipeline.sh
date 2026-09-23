#!/bin/bash
#SBATCH --job-name=inter
#SBATCH --output=%x.log
#SBATCH --error=%x.error
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G


# Correr pipeline de nextflow para encontrar interacciones, (metodo binomial, todos los sitios) 

date
echo "===== Beginning pipeline ====="

# cargar modulos
module load nextflow/main
module load R/4.5.3

# comands here
nextflow -version
nextflow run MAGs_BGCs_interactions/pipeline.nf \
  --indir metadata/ \
  --outdir interactions/ \
  --microbial_lineage mOTUs_Species_Cluster \
  --bgc_groups gcc \
  --method binomial \
  --quality 1

echo "===== Pipeline done ====="
date
