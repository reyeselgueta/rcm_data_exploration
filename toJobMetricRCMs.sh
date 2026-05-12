#!/bin/bash
#SBATCH --job-name=metric-process
#SBATCH --time=14:00:00
#SBATCH --mem=16Gb
#SBATCH --partition=wncompute_meteo
#SBATCH --output=outs/metric_process_%j.out
#SBATCH --error=outs/metric_process_%j.err

source ~/.bashrc
conda activate deep4downscaling-gpu
cd /gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration

target_coords=$1
frecuency=$2
statistic=$3
scenario=$4
variable=$5
spatial=$6

python metric_process_rcms.py $target_coords $frecuency $statistic $scenario $variable $spatial
# sbatch toJobMetricProcess.sh valencia 1hr mean historical pr True
