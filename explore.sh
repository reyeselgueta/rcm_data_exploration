#!/bin/bash
#SBATCH --job-name=explore
#SBATCH --output=outs/explore-%j.out
#SBATCH --error=outs/explore-%j.err
#SBATCH --partition=meteo_long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60Gb
#SBATCH --time=12:00:00
#SBATCH --mail-user=reyes@ifca.unican.es
#SBATCH --mail-type=END,FAIL


source ~/miniconda3/etc/profile.d/conda.sh
conda activate cdo-simple
cd /nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration


python explore.py

