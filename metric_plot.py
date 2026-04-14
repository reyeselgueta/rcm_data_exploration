#!git clone https://github.com/SantanderMetGroup/deep4downscaling.git
#sbatch --export target_coords='valencia',period=2,frecuency='day' toJobPreprocess.sh
import time
import xarray as xr
import numpy as np
from pathlib import Path
import sys
import utils_explore as utils
from utils_explore import rcm_dict, frecuency_dict, value_limits, colobar_limits, mask_dict
import os
import glob
from utils_dana import TimeSeriesPlot, colors

BASE_PATH = Path("/gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration/deep4downscaling")#"/vols/abedul/home/meteo/reyess/paper1-code/deep4downscaling")
sys.path.insert(0, str(BASE_PATH))

import deep4downscaling.viz
import deep4downscaling.trans
import deep4downscaling.metrics
import deep4downscaling.metrics_ccs


DATA_PATH = './data/input'
FIGURES_PATH = '/gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration/example_figures/'
MODELS_PATH = './models'
ASYM_PATH = './data/asym'
DATA_PATH_METRICS = '/gpfs/projects/meteo/WORK/reyesjsf/data/rcm-gcm/'

# target_coords = sys.argv[1] # 'valencia' 'iberia' 'europe'
# statistic = sys.argv[2] # 'max' or 'mean' or 'date_max''
# spatial = sys.argv[3] # 'spatial' or 'aggregated'
# per_grade = sys.argv[4] # 'True' or 'False'

target_coords = 'valencia'
statistic = 'mean'
spatial = 'False'
per_grade = 'False'


spatial = True if spatial == 'True' else False
per_grade = True if per_grade == 'True' else False
lat_center, lon_center, folder_name, grid_multiplier = utils.get_coords(target_coords)
years_hist = ('1994', '2014') 
years_future = ('2080', '2100')
seasons = ['Annual', 'DJF', 'MAM', 'JJA', 'SON']
frecuencies = ['day', '1hr']
rcms = ['BCCR-UCAN', 'BCCR-UCAN_eur12', 'CNRM-MF']
time_start = time.time()

if spatial == True:
    # Dictionary to store metrics for plot
    data_to_plot = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in rcms} for season in seasons}

    # We go over RCMs to load data
    for frecuency in frecuencies:
        for rcm_name in rcms:
            # LOad data
            target_path = f"{DATA_PATH_METRICS}{folder_name}/climatology_pr_{frecuency}_{rcm_name}_{target_coords}_{statistic}_{years_future[0]}-{years_future[1]}.nc"
            secundary_path = f"{DATA_PATH_METRICS}{folder_name}/climatology_pr_{frecuency}_{rcm_name}_{target_coords}_{statistic}_{years_hist[0]}-{years_hist[1]}.nc"
            ds_target = xr.open_dataset(target_path)
            ds_secundary = xr.open_dataset(secundary_path)
            if per_grade == True and rcm_name != 'CNRM-MF':
                tas_fut_path = f"{DATA_PATH_METRICS}{folder_name}/climatology_tas_{frecuency}_{rcm_name}_{target_coords}_mean_{years_future[0]}-{years_future[1]}.nc"
                tas_hist_path = f"{DATA_PATH_METRICS}{folder_name}/climatology_tas_{frecuency}_{rcm_name}_{target_coords}_mean_{years_hist[0]}-{years_hist[1]}.nc"
                ds_tas_fut = xr.open_dataset(tas_fut_path)
                ds_tas_hist = xr.open_dataset(tas_hist_path)
            
            for season in seasons:
                total_percentage = (ds_target['pr'].sel(season=season) - ds_secundary['pr'].sel(season=season)) / ds_secundary['pr'].sel(season=season) * 100
                if per_grade == True:
                    tas_diff = ds_tas_fut['tas'].mean(dim='time') - ds_tas_hist['tas'].mean(dim='time') if rcm_name != 'CNRM-MF' else 3
                    total_percentage = total_percentage / tas_diff
                data_to_plot[season][rcm_name][frecuency] = total_percentage

    for season in seasons:
        utils.multi_map(data=data_to_plot[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['relative_per_degree' if per_grade else 'relative'],
                    color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['relative_per_degree' if per_grade else 'relative'],
                    title=f'Relative (per degree) {statistic}-{season} CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12 {'per °' if per_grade else ''}',
                    fig_path=FIGURES_PATH, fig_name=f'Relative_{statistic}_{season}_{target_coords}{'_per °' if per_grade else ''}.png')
                    #mask_dict = mask_dict)
else:
    print("Generating timeseries plots for spatially aggregated metrics")
    data_to_plot = {frecuency: [] for frecuency in frecuencies}
    for frecuency in frecuencies:
        # data_to_plot = {frecuency: {rcm_name: None for rcm_name in rcms} for frecuency in frecuencies}
        for rcm_name in rcms:
            target_path = f"{DATA_PATH_METRICS}{folder_name}/climatology_pr_{frecuency}_{rcm_name}_{target_coords}_{statistic}{'_spatial' if spatial==True else ''}_{years_future[0]}-{years_future[1]}.nc"
            data_to_plot[frecuency].append({'data_list':[xr.open_dataset(target_path).pr], 'label': rcm_name, 'color': colors[rcm_name]})


    Plot = TimeSeriesPlot(n_rows=len(frecuencies), n_cols=1, figsize_scale=4)

    for var_idx, frecuency in enumerate(frecuencies):
        time_var = 'time'
        # ===========================
        # PLOT
        # ===========================
        Plot.plot_subplot(
            row_idx=var_idx,
            col_idx=0,
            data_metric=data_to_plot[frecuency],
            time_var=time_var,
            title=f"Algo",
            y_label=f"mm/{frecuency}",
            x_label=frecuency.capitalize()
        )

    Plot.finalize(savepath=f'{FIGURES_PATH}{folder_name}/timeseries_pr_{target_coords}_{statistic}_{years_future[0]}-{years_future[1]}.png', title=f'Timeseries of metrics for {frecuency} {statistic} {target_coords}')
    Plot.finalize(savepath=f'{FIGURES_PATH}{folder_name}/timeseries_pr_{target_coords}_{statistic}_{years_future[0]}-{years_future[1]}.pdf', title=f'Timeseries of metrics for {frecuency} {statistic} {target_coords}')


time_end = time.time()
print(f"Plots generados para {frecuency} : {(time_end-time_start)/60:.2f} minutes")


