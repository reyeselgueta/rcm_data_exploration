#!/bin/bash
#SBATCH --job-name=metric-plot
#SBATCH --time=14:00:00
#SBATCH --mem=8Gb
#SBATCH --partition=wncompute_meteo
#SBATCH --output=outs/metric_plot_%j.out
#SBATCH --error=outs/metric_plot_%j.err

source ~/.bashrc
conda activate deep4downscaling-gpu
cd /gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration

target_coords=$1
statistic=$2
spatial=$3
per_grade=$4

python metric_plot.py $target_coords $statistic $spatial $per_grade
# sbatch toJobMetricProcess.sh valencia 1hr mean historical pr True
