#!git clone https://github.com/SantanderMetGroup/deep4downscaling.git
#sbatch --export target_coords='valencia',period=2,frecuency='day' toJobPreprocess.sh
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
period = int(sys.argv[2]) # 1 or 2
frecuency = sys.argv[3] #1hr or day
statistic = sys.argv[4] # 'max' or 'mean' or 'date_max''
# target_coords = 'valencia' 
# period = 1
# statistic = 'date_max'
# frecuency = 'day'
print(f"Target coords: {target_coords} Period: {period} Frecuency: {frecuency} Statistic: {statistic}")



# EURR-3 solo tiene historical y de ICTP

if target_coords == 'iberia':
    lat_center = 39.000
    lon_center = -2.875
    folder_name = 'iberia-3'
    grid_multiplier = 3
elif target_coords == 'valencia':
    lat_center = 39.467
    lon_center = -0.375
    folder_name = 'val-3'
    grid_multiplier = 1
elif target_coords == 'europe':
    lat_center = 39.467
    lon_center = -0.375
    folder_name = 'eur-3'
    grid_multiplier = 9
else:
    raise ValueError("target_coords must be 'iberia', 'valencia', or 'europe'")

seasons = ['Annual', 'DJF', 'MAM', 'JJA', 'SON']

time1= time.time()
for rcm_name, data_rcm in utils.data_rcm_paths.items():
    # General parameters
    grid_number = 90 if rcm_name != 'BCCR-UCAN_eur12' else 30
    grid_number = grid_number * grid_multiplier
    gcm_name = data_rcm['gcm']
    time2 = time.time()
    time_multiplier = 3600 if frecuency == '1hr' else 86400
    #FOR CLIMATOLOGY
    if period == 1:
        data_paths = data_rcm['historical']
        utils.mem_status("Antes de abrir datos hist")
        data = xr.open_mfdataset(
                f"{data_paths[frecuency]}*.nc", 
                combine="by_coords",
                chunks={'time': 100}
            ).sel(time=slice('1995', '2014'))

        lat=data.lat.values
        lon=data.lon.values
        abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
        iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

        # Slice area around the city
        half_box = grid_number // 2
        y_slice = slice(max(0, iy - half_box), iy + half_box)
        x_slice = slice(max(0, ix - half_box), ix + half_box)

        data_selected = data.pr.isel(y=y_slice, x=x_slice)
        data_climatology = {season: None for season in seasons}

        for season in seasons:
            data_yearly = utils.get_season_data(data_selected, season, year_to_drop=2015)
            data_statistic = data_statistic = utils.get_statistic(data_yearly, statistic, time_multiplier, season)
            data_climatology[season] = data_statistic# - data_prh_climatology[season][rcm_name][gcm_name]

        ds_combined = xr.concat(
            [data_climatology[season] for season in seasons],
                dim='season'
            )

        ds_combined = ds_combined.assign_coords(season=seasons)
        ds_combined.to_netcdf(f"{DATA_PATH_METRICS}{folder_name}/climatology_{frecuency}_{rcm_name}_{target_coords}_{statistic}_1995-2015.nc", compute=True)


        time3 = time.time()
        print(f"Datos guardados hist {frecuency} {rcm_name}: {(time3-time2)/60:.2f} minutes")

    # SSP370
    if period == 2:
        if data_rcm['ssp370'] is None:
            continue
        data_paths = data_rcm['ssp370'][frecuency]
        years_fut = utils.gcm_gwl3_years[gcm_name]
        files = glob.glob(f"{data_paths}*.nc")
        files_selected = []
        years = range(int(years_fut[0]), int(years_fut[1])+1)
        for f in files:
            if any(str(year) in f for year in years):
                files_selected.append(f)
        if not files_selected:
            raise ValueError(f"No files for Rcm: {rcm_name} Frecuency:{frecuency} Period:{period} Years: {years_fut}")
        # Abrir solo los que interesan
        data = xr.open_mfdataset(files_selected,
                        combine="by_coords",
                        chunks={'time': 100}
                    ).sel(time=slice(years_fut[0], years_fut[1]))

        lat=data.lat.values
        lon=data.lon.values
        abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
        iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

        # Slice area around the city
        half_box = grid_number // 2
        y_slice = slice(max(0, iy - half_box), iy + half_box)
        x_slice = slice(max(0, ix - half_box), ix + half_box)

        data_selected = data.pr.isel(y=y_slice, x=x_slice)
        data_future = {season: {rcm_name : {gcm_name: gcm_name}} for season in seasons}

        for season in seasons:
            data_yearly = utils.get_season_data(data_selected, season, year_to_drop=int(years_fut[1])+1)
            data_statistic = utils.get_statistic(data_yearly, statistic, time_multiplier, season)
            data_future[season] = data_statistic #- data_climatology[season][rcm_name][data_rcm['gcm']]

        ds_combined = xr.concat(
            [data_future[season] for season in seasons],
                dim='season'
            )
        ds_combined = ds_combined.assign_coords(season=seasons)
        ds_combined.to_netcdf(f"{DATA_PATH_METRICS}{folder_name}/ssp370_{frecuency}_{rcm_name}_{target_coords}_{statistic}_CWL3.nc", compute=True)

        time3= time.time()
        print(f"Datos guardados future {frecuency} {rcm_name}: {(time3-time2)/60:.2f} minutes")



time5 = time.time()
print(f"Datos guardados: {(time5-time1)/60:.2f} minutes")