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



ERA5_PATH = '/lustre/gmeteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/download/ERA5/pr/*.nc'
ERA5_LAND_PATH = '/lustre/gmeteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/download/ERA5-land/pr/*.nc'
CERRA_PATH = '/gpfs/projects/meteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/CERRA/final_products/climate_index/CERRA/pr/gr006/'
CERRA_LAND_PATH = '/gpfs/projects/meteo/WORK/PROYECTOS/2022_C3S_Atlas/workflow/datasets/CICAv2/CERRA-land/final_products/Global/CERRA-land/year/cdd/gr006/day/'
EOBS_PATH = '/lustre/gmeteo/WORK/DATA/C3S-CDS/CDS-Curated-Data/raw/insitu-gridded-observations-europe/daily/native/rr/rr_insitu-gridded-observations-europe_1950-2024_31_0e.nc'
ROCIO_PATH = '/gpfs/projects/meteo/WORK/reyesjsf/data/rocio/'
