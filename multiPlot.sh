#!/bin/bash


target_coords_list=(
    "valencia"
    #"iberia"
    #"europe"
)



for target_coords in "${target_coords_list[@]}"; do
    for custom in {1..1}; do
        echo "Submitting job for Coords: $target_coords custom: $custom"
        sbatch --export=ALL,target_coords=$target_coords,plot_num='5',flux=0,custom=$custom toJobPlot.sh
    done
done