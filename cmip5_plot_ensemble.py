import time
import xarray as xr
import numpy as np
from pathlib import Path
import sys
import utils_explore as utils


DATA_PATH_METRICS = '/nfs/home/gmeteo/reyess/rcm_exploration/data/metrics/'
FIG_PATH = '/nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration/example_figures/relatives/'
# metric = sys.argv[1] # mon-mean, yr-mean, max-mean
metric = "yr-mean" # mon-mean, yr-mean, max-mean
area = "Valencia"

lat_target, lon_target = utils.get_slice_coords(area)


rcm_list = ['CCLM4-8-17', 'COSMO-crCLIM-v1-1', 'HIRHAM5',
            'HadREM3-GA7-05', 'RACMO22E']
gcm_list = ['CNRM-CM5', 'EC-EARTH', 'HadGEM2-ES', 'MPI-ESM-LR', 'NorESM1-M']
seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']
row_metrics = ['Mean', 'P80', 'P20']


rcm_dict = {rcm_name:rcm_name for rcm_name in rcm_list}
gcm_dict = {gcm_name:gcm_name for gcm_name in gcm_list}

gcm_gwl3_years = utils.gcm_gwl3_years
prhmax_data = {row_metric: {season: None for season in seasons} for row_metric in row_metrics}
rx1day_data = {row_metric: {season: None for season in seasons} for row_metric in row_metrics}

# Load metrics
for row_metric in row_metrics:
    ds_rx1day = xr.open_dataset(f'{DATA_PATH_METRICS}/cmip5ensemble_relative_{row_metric.lower()}_{metric}_rx1day.nc')
    ds_prhmax = xr.open_dataset(f'{DATA_PATH_METRICS}/cmip5ensemble_relative_{row_metric.lower()}_{metric}_prhmax.nc')

    for season in seasons:
        rx1day_data[row_metric][season] = ds_rx1day.sel(season=season)
        prhmax_data[row_metric][season] = ds_prhmax.sel(season=season)


utils.multi_map(data=rx1day_data, x_map=row_metrics, y_map=seasons, vlimits=[(0, 30), (0, 30), (0, 30)], var='rx1day',
        color=['BrBG', 'BrBG', 'BrBG'], cbar_limits=[(0, 10, 10), (0, 10, 10), (0, 10, 10)], title=f'Metrics rx1day {metric}',
        fig_path=FIG_PATH, fig_name=f'Metric_Rx1day_{metric}.png')

utils.multi_map(data=prhmax_data, x_map=row_metrics, y_map=seasons, vlimits=[(0, 30), (0, 30), (0, 30)], var='prhmax',
        color=['BrBG', 'BrBG', 'BrBG'], cbar_limits=[(0, 10, 10), (0, 10, 10), (0, 10, 10)], title=f'Metrics prhmax {metric}',
        fig_path=FIG_PATH, fig_name=f'Metric_prhmax_{metric}.png')