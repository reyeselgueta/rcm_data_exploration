#!/bin/bash
#SBATCH --job-name=data-process
#SBATCH --time=14:00:00
#SBATCH --mem=48Gb
#SBATCH --partition=wncompute_meteo
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

source ~/.bashrc
conda activate deep4downscaling-gpu
cd /gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration
python data_process.py $target_coords $period $frecuency $statistic
