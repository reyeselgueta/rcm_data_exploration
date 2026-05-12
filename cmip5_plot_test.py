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

lat_iberia = slice(34.12, 43.88)
lon_iberia = slice(-9.875, 5.875)
lat_valencia = slice(38, 41)
lon_valencia = slice(-2.5, 2.0)


rcm_list = ['CCLM4-8-17', 'COSMO-crCLIM-v1-1', 'HIRHAM5',
            'HadREM3-GA7-05', 'RACMO22E']
gcm_list = ['CNRM-CM5', 'EC-EARTH', 'HadGEM2-ES', 'MPI-ESM-LR', 'NorESM1-M']

rcm_dict = {rcm_name:rcm_name for rcm_name in rcm_list}
gcm_dict = {gcm_name:gcm_name for gcm_name in gcm_list}

gcm_gwl3_years = utils.gcm_gwl3_years
relative_diff = {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list}

for rcm_name in rcm_list:
    for gcm_name in gcm_list:
        # Load metrics
        path_rx1day_mean = f'{DATA_PATH_METRICS}/cmip5ensemble_relative_mean_{metric}_rx1day.nc'
        path_rx1day_p20 = f'{DATA_PATH_METRICS}/cmip5ensemble_relative_p20_{metric}_rx1day.nc'
        path_rx1day_p80 = f'{DATA_PATH_METRICS}/cmip5ensemble_relative_p80_{metric}_rx1day.nc'

        path_prhmax_mean = f'{DATA_PATH_METRICS}/cmip5ensemble_relative_mean_{metric}_prhmax.nc'
        path_prhmax_p20 = f'{DATA_PATH_METRICS}/cmip5ensemble_relative_p20_{metric}_prhmax.nc'
        path_prhmax_p80 = f'{DATA_PATH_METRICS}/cmip5ensemble_relative_p80_{metric}_prhmax.nc'
utils.multi_map(data=relative_diff, x_map=rcm_dict, y_map=gcm_dict, vlimits=(-50, 50), var='prhmax',
        color='BrBG', cbar_limits=(0, 10, 10), title=f'Relative difference prhmax - rx1day ({metric}) - % (Average rx1day over 20 years, 1986-2005 as reference, and gwl3 as target.))',
        fig_path=FIG_PATH, fig_name=f'Difference_Relatives_{metric}.png')