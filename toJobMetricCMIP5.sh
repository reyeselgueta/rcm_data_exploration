#!/bin/bash
#SBATCH --job-name=cmip5-metrics
#SBATCH --output=outs/cmip5-metrics-%j.out
#SBATCH --error=outs/cmip5-metrics-%j.err
#SBATCH --partition=meteo_long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=wn055
#SBATCH --mem=40Gb
#SBATCH --time=12:00:00
#SBATCH --mail-user=reyes@ifca.unican.es
#SBATCH --mail-type=END,FAIL


source ~/miniconda3/etc/profile.d/conda.sh
conda activate cdo-simple
cd /nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration

target_var=$1
metric=$2

python cmip5_metrics.py $target_var $metric
# sbatch toJobMetricProcess.sh valencia 1hr mean historical pr True
