#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --time=14:00:00
#SBATCH --mem=16Gb
#SBATCH --partition=wncompute_meteo
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

source ~/.bashrc
conda activate deep4downscaling-gpu
cd /gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration
python plot.py $target_coords $plot_num $flux $custom
