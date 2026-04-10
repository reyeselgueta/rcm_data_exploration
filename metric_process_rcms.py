#!git clone https://github.com/SantanderMetGroup/deep4downscaling.git
#sbatch --export target_coords='valencia',period=2,frecuency='day' toJobPreprocess.sh

# This script has the objective of processing the RCM data to obtain the climatology of the selected statistic (max, mean, date_max, P95, P995) for each season and save it in a netcdf file. This is done for each RCM and GCM combination, and for the selected target coordinates (valencia, iberia, europe). The processed data is saved in the specified path for metrics, with a filename that includes the target coordinates, statistic, frequency, and years of the scenario. The script uses xarray to handle the netcdf files and perform the necessary calculations. The time taken for processing is also printed at the end.
# And also for percetange of 
import time
import xarray as xr
import numpy as np
from pathlib import Path
import sys
import utils_explore as utils
import os
import glob

BASE_PATH = Path("/gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration/deep4downscaling")#"/vols/abedul/home/meteo/reyess/paper1-code/deep4downscaling")
sys.path.insert(0, str(BASE_PATH))

import deep4downscaling.viz
import deep4downscaling.trans
import deep4downscaling.metrics
import deep4downscaling.metrics_ccs


DATA_PATH = './data/input'
FIGURES_PATH = '/gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration/deep4downscaling/notebooks/figures'
MODELS_PATH = './models'
ASYM_PATH = './data/asym'
DATA_PATH_METRICS = '/gpfs/projects/meteo/WORK/reyesjsf/data/rcm-gcm/'


target_coords = sys.argv[1] # 'valencia' 'iberia' 'europe'
frecuency = sys.argv[2] #1hr or day
statistic = sys.argv[3] # 'max' or 'mean' or 'date_max''
scenario = sys.argv[4] # 'historical' or 'ssp370'
variable = sys.argv[5] # 'pr' or 'tas'
spatial = sys.argv[6] # 'spatial' or 'aggregated'
# target_coords = 'valencia'
# frecuency = 'day'
# statistic = 'mean-max'
# scenario = 'historical' or 'ssp370'
# variable = sys.argv[5] # 'pr' or 'tas'
# spatial = sys.argv[6] # 'spatial' or 'aggregated'

spatial = True if spatial == 'True' else False


time_start = time.time()
print(f"Target coords: {target_coords} Frecuency: {frecuency} Statistic: {statistic} Variable: {variable} Scenario: {scenario} Spatial: {spatial}")
seasons = ['Annual', 'DJF', 'MAM', 'JJA', 'SON']

lat_center, lon_center, folder_name, grid_multiplier = utils.get_coords(target_coords)

# Spatial metrics:

for rcm_name, data_rcm in utils.data_rcm_paths.items():
    # General parameters
    grid_number = 90 if rcm_name != 'BCCR-UCAN_eur12' else 30
    grid_number = grid_number * grid_multiplier
    gcm_name = data_rcm['gcm']
    
    time_multiplier = 3600 if frecuency == '1hr' else 86400
    time_multiplier = 1 if variable != 'pr' else time_multiplier


    if data_rcm[scenario] is None:
        print(f"No data for Rcm: {rcm_name} Scenario: {scenario}")
        continue
    if data_rcm[scenario][variable] is None:
        print(f"No data for Rcm: {rcm_name} Variable: {variable}")
        continue
    data_paths = data_rcm[scenario][variable][frecuency]
    #years_scenario = utils.gcm_gwl3_years[gcm_name]
    years_scenario = ('1994', '2014') if scenario == 'historical' else ('2080', '2100')
    files = glob.glob(f"{data_paths}*.nc")
    files_selected = []
    years = range(int(years_scenario[0]), int(years_scenario[1])+1)
    for f in files:
        if any(str(year) in f for year in years):
            files_selected.append(f)
    if not files_selected:
        raise ValueError(f"No files for Rcm: {rcm_name} Frecuency:{frecuency} Years: {years_scenario} Scenario: {scenario} Spatial: {spatial} Variable {variable}")
    # Open only the selected files and concatenate them along the time dimension, selecting only the years of interest
    data = xr.open_mfdataset(files_selected,
                    combine="by_coords",
                    chunks={'time': 100}
                ).sel(time=slice(years_scenario[0], years_scenario[1]))
    print(f"Data loaded for Rcm: {rcm_name} Frecuency:{frecuency} Years: {years_scenario} Scenario: {scenario} Spatial: {spatial} Variable {variable}")
    print(data)

    # Find the closest grid point to the target coordinates
    lat=data.lat.values
    lon=data.lon.values
    abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
    iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

    # Slice area around the city
    half_box = grid_number // 2
    y_slice = slice(max(0, iy - half_box), iy + half_box)
    x_slice = slice(max(0, ix - half_box), ix + half_box)

    data_selected = data[variable].isel(y=y_slice, x=x_slice)
    if spatial == True:
        data_climatology =  {season: None for season in seasons} #{season: {rcm_name : {gcm_name: gcm_name}} for season in seasons}

        for season in seasons:
            data_yearly = utils.get_season_data(data_selected, season, year_to_drop=2015)
            print(f"Data yearly in season: {season}")
            print(data_yearly)
            data_statistic = utils.get_statistic(data_yearly, statistic, time_multiplier, season, spatial=spatial)
            data_climatology[season] = data_statistic# - data_prh_climatology[season][rcm_name][gcm_name]

        ds_combined = xr.concat(
            [data_climatology[season] for season in seasons],
            dim="season")#xr.DataArray(seasons, dims="season", name="season")
        ds_combined = ds_combined.assign_coords(season=seasons)
    elif spatial == False:
        data_yearly = utils.get_season_data(data_selected, 'Annual', year_to_drop=2015)
        ds_combined = utils.get_statistic(data_yearly, statistic, time_multiplier, 'Annual', spatial=spatial)

    ds_combined.to_netcdf(f"{DATA_PATH_METRICS}{folder_name}/climatology_{variable}_{frecuency}_{rcm_name}_{target_coords}_{statistic}_{years_scenario[0]}-{years_scenario[1]}.nc", compute=True)


time_end = time.time()
print(f"Data saved for {variable} {frecuency} {target_coords} {statistic} {scenario}: {(time_end-time_start)/60:.2f} minutes")