import time
import xarray as xr
import cf_xarray as cfxr
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



FIGURES_PATH = '/gpfs/users/reyesjsf/rcm-exploration/rcm_data_exploration/deep4downscaling/notebooks/figures'
DATA_PATH_METRICS = '/gpfs/projects/meteo/WORK/reyesjsf/data/rcm-gcm/'
DATA_PATH_METRICS = '/gpfs/projects/meteo/WORK/reyesjsf/data/rcm-gcm/'



ERA5_PATH = '/lustre/gmeteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/download/ERA5/pr/*.nc'
ERA5_LAND_PATH = '/lustre/gmeteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/download/ERA5-land/pr/*.nc'
CERRA_PATH = '/gpfs/projects/meteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/CERRA/download/CERRA/pr/*.nc'
CERRA_LAND_PATH = '/gpfs/projects/meteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/CERRA-land/download/Global/CERRA-Land/pr/day/*.nc'
EOBS_PATH = '/lustre/gmeteo/WORK/DATA/C3S-CDS/CDS-Curated-Data/raw/insitu-gridded-observations-europe/daily/native/rr/rr_insitu-gridded-observations-europe_1950-2024_31_0e.nc'
ROCIO_PATH = '/lustre/gmeteo/WORK/reyess/data/predictand/AEMET_0.25deg/pr/*.nc'


lon = (-1.75, 1.05)
lat = (38.45, 40.55)
lon_min, lon_max = lon[0], lon[1]
lat_min, lat_max = lat[0], lat[1]

# CERRA Dana day
ds_cerra = xr.open_mfdataset(CERRA_PATH, combine='by_coords')
ds_cerra_dana = ds_cerra.sel(valid_time="2024-10-29")

mask = (
    (ds_cerra_dana.longitude >= lon_min) &
    (ds_cerra_dana.longitude <= lon_max) &
    (ds_cerra_dana.latitude >= lat_min) &
    (ds_cerra_dana.latitude <= lat_max)
)
mask = mask.compute()  # Compute the mask to avoid lazy evaluation issues

ds_cerra_dana = ds_cerra_dana.where(mask, drop=True)

# CERRA LAND Dana day
ds_cerra_land = xr.open_mfdataset(CERRA_LAND_PATH, combine='by_coords')
ds_cerra_land_dana = ds_cerra_land.sel(valid_time="2024-10-29T06:00:00.000000000", method='nearest')

mask = (
    (ds_cerra_land_dana.longitude >= lon_min) &
    (ds_cerra_land_dana.longitude <= lon_max) &
    (ds_cerra_land_dana.latitude >= lat_min) &
    (ds_cerra_land_dana.latitude <= lat_max)
)
mask = mask.compute()  # Compute the mask to avoid lazy evaluation issues

ds_cerra_land_dana = ds_cerra_land_dana.where(mask, drop=True)



