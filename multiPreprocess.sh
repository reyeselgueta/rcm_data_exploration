#!/bin/bash


target_coords_list=(
    "valencia"
    #"iberia"
    #"europe"
)

frecuency_list=(
    "day"
    "1hr"
)

statistic_list=(
    #"mean"
    #"max"
    "date_max"
)

for target_coords in "${target_coords_list[@]}"; do
    for period in {1..2}; do
        for frecuency in "${frecuency_list[@]}"; do
            for statistic in "${statistic_list[@]}"; do
                echo "Submitting job for Coords: $target_coords, Period: $period, Frecuency: $frecuency, Statistic: $statistic"
                sbatch --export=ALL,target_coords=$target_coords,period=$period,frecuency=$frecuency,statistic=$statistic toJobPreprocess.sh
            done
        done
    done
done