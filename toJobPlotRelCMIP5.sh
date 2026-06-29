#!/bin/bash
#SBATCH --job-name=cmip5-plot
#SBATCH --output=outs/cmip5-plot-%j.out
#SBATCH --error=outs/cmip5-plot-%j.err
#SBATCH --partition=meteo_long
#SBATCH --nodes=1
#SBATCH --exclude=wn055
#SBATCH --ntasks=1
#SBATCH --mem=10Gb
#SBATCH --time=12:00:00
#SBATCH --mail-user=reyes@ifca.unican.es
#SBATCH --mail-type=END,FAIL

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cdo-simple
cd /nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration

metric=$1
area=$2

python cmip5_plot_relative.py $metric $area

