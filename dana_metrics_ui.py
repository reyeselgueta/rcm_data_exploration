import time
import xarray as xr
import numpy as np
from pathlib import Path
import sys
import utils_explore as utils
import os
import glob

# BASE_PATH = Path("/gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration/deep4downscaling")#"/vols/abedul/home/meteo/reyess/paper1-code/deep4downscaling")
# sys.path.insert(0, str(BASE_PATH))

# import deep4downscaling.viz
# import deep4downscaling.trans
# import deep4downscaling.metrics
# import deep4downscaling.metrics_ccs


DATA_PATH = './data/input'
FIGURES_PATH = '/nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration/example_figures/'
MODELS_PATH = './models'
ASYM_PATH = './data/asym'
DATA_PATH_METRICS = '/nfs/home/gmeteo/reyess/rcm_exploration/data/metrics/'



ERA5_PATH = '/lustre/gmeteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/download/ERA5/pr/'
ERA5_LAND_PATH = '/lustre/gmeteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/download/ERA5-land/pr/'
CERRA_PATH = '/gpfs/projects/meteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/CERRA/download/CERRA/pr/*.nc'
CERRA_LAND_PATH = '/gpfs/projects/meteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/CERRA-land/download/Global/CERRA-Land/pr/day/*.nc'
EOBS_PATH = '/lustre/gmeteo/WORK/DATA/C3S-CDS/CDS-Curated-Data/raw/insitu-gridded-observations-europe/daily/native/rr/rr_insitu-gridded-observations-europe_1950-2024_31_0e.nc'
ROCIO_PATH = '/lustre/gmeteo/WORK/reyess/data/predictand/AEMET_0.25deg/pr/'


lon = (-1.75, 1.05)
lat = (38.45, 40.55)
lon_min, lon_max = lon[0], lon[1]
lat_min, lat_max = lat[0], lat[1]

# ERA5 Dana day (Dato horario)
ds_era5 = xr.open_dataset(f'{ERA5_PATH}total_precipitation-reanalysis-2024-01-01_2024-12-31.nc')
#ds_era5 = xr.open_mfdataset(f"ERA5_PATH}*.nc", combine="nested", concat_dim="time")#combine='by_coords')
ds_era5_dana = ds_era5.sel(valid_time="2024-10-29")

mask = (
    (ds_era5_dana.longitude >= lon_min) &
    (ds_era5_dana.longitude <= lon_max) &
    (ds_era5_dana.latitude >= lat_min) &
    (ds_era5_dana.latitude <= lat_max)
)
mask = mask.compute()  # Compute the mask to avoid lazy evaluation issues

ds_era5_dana = ds_era5_dana.where(mask, drop=True)
ds_era5_dana_mean = ds_era5_dana.mean(dim='valid_time')
ds_era5_dana_mean.to_netcdf(f"{DATA_PATH_METRICS}/ERA5_dana_mean_hourly_2024-10-29.nc")


# ERA5 LAND Dana day
#ds_era5_land = xr.open_mfdataset(f"{ERA5_LAND_PATH}*.nc", combine='by_coords')
ds_era5_land = xr.open_dataset(f'{ERA5_LAND_PATH}total_precipitation-2024-10-01_2024-10-31.nc')
ds_era5_land_dana = ds_era5_land.sel(valid_time="2024-10-29", method='nearest')

mask = (
    (ds_era5_land_dana.longitude >= lon_min) &
    (ds_era5_land_dana.longitude <= lon_max) &
    (ds_era5_land_dana.latitude >= lat_min) &
    (ds_era5_land_dana.latitude <= lat_max)
)
mask = mask.compute()  # Compute the mask to avoid lazy evaluation issues

ds_era5_land_dana = ds_era5_land_dana.where(mask, drop=True)
ds_era5_land_dana.to_netcdf(f"{DATA_PATH_METRICS}/ERA5-Land_dana_mean_daily_2024-10-29.nc")

# EOBS Dana day
ds_eobs = xr.open_dataset(f'{EOBS_PATH}')
ds_eobs_dana = ds_eobs.sel(time="2024-10-29")

mask = (
    (ds_eobs_dana.longitude >= lon_min) &
    (ds_eobs_dana.longitude <= lon_max) &
    (ds_eobs_dana.latitude >= lat_min) &
    (ds_eobs_dana.latitude <= lat_max)
)
mask = mask.compute()  # Compute the mask to avoid lazy evaluation issues

ds_eobs_dana = ds_eobs_dana.where(mask, drop=True)
ds_eobs_dana.to_netcdf(f"{DATA_PATH_METRICS}/EOBS_dana_mean_daily_2024-10-29.nc")


# ROCIO Dana day
ds_rocio = xr.open_mfdataset(ROCIO_PATH, combine='by_coords')
ds_rocio_dana = ds_rocio.sel(time="2024-10-29")


