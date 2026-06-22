#!/bin/bash

#SBATCH --job-name=6_go_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6_go_enrichment.err"
#SBATCH --output="./logs/6_go_enrichment.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-gsea

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"

cd $BASE_DIR

# Run the R script
Rscript scripts/6_go_enrichment.R
