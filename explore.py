!git clone https://github.com/SantanderMetGroup/deep4downscaling.git


import xarray as xr
import numpy as np
from pathlib import Path
import sys
import utils_explore as utils

BASE_PATH = Path("...")
sys.path.insert(0, str(BASE_PATH))

import deep4downscaling.viz
import deep4downscaling.trans
import deep4downscaling.metrics
import deep4downscaling.metrics_ccs

DATA_PATH = './data/input'
FIGURES_PATH = '/notebooks/figures'
MODELS_PATH = './models'
ASYM_PATH = './data/asym'
DATA_PATH_PREDICTAND = '....'

data_url = '...'
data_fut_url = '...'
target_var = 'rx1day'
data_prh_url = '...'


lat_iberia = slice(34.12, 43.88)
lon_iberia = slice(-9.875, 5.875)
lat_valencia = slice(38, 41)
lon_valencia = slice(-2.5, 2.0)


# Create general structure
rcm_list = ['ALADIN63', 'CCLM4-8-17', 'COSMO-crCLIM-v1-1', 'HIRHAM5',
            'HadREM3-GA7-05', 'RACMO22E', 'RCA4', 'REMO2009',
            'REMO2015', 'RegCM4-6', 'WTF361H', 'WRF381P']
gcm_list = ['CNRM-CM5', 'CanESM2', 'EC-EARTH', 'HadGEM2-ES',
            'IPSL-CM5A-MR', 'MIROC5', 'MPI-ESM-LR', 'NorESM1-M']

rcm_dict = {rcm_name:rcm_name for rcm_name in rcm_list}
gcm_dict = {gcm_name:gcm_name for gcm_name in gcm_list}

# CLIMATOLOGY DX1DAY
data_hist = xr.open_dataset(f'{data_url}historical/{target_var}_CORDEX-EUR-11_historical_mon_197001-200512_v02.nc')
data_hist_ref = data_hist.sel(time=slice('1986','2005'))
data_hist_ref = data_hist_ref.sel(lat=lat_iberia, lon=lon_iberia)#IBERIA
data_hist_ref_r1 = data_hist_ref.where(data_hist_ref.gcm_variant == 'r1i1p1', drop=True)
data_dx_climatology = {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list}
for rcm_name in rcm_list:
    data_hist_rcm = data_hist_ref_r1.where(data_hist_ref_r1.rcm_model==rcm_name, drop=True)
    for gcm_name in gcm_list:
        data_hist_selected = data_hist_rcm.where(data_hist_rcm.gcm_model==gcm_name, drop=True)
        if data_hist_selected.sizes.get("member", 0) == 0:
            continue
        data_hist_selected_mean = data_hist_selected.mean(dim='time').squeeze('member')
        data_dx_climatology[rcm_name][gcm_name] = data_hist_selected_mean

utils.multi_map(data=data_dx_climatology, x_map=rcm_dict, y_map=gcm_dict, vlimits=(0, 80), var='rx1day',
        color='Blues', cbar_limits=(0, 10, 10), title='Climatology rx1day - mm (Average rx1day over 20 years, 1986-2005)',
        fig_path=FIGURES_PATH, fig_name='climatology_rx1day.png')

# DELTA DX1DAY
data_fut = xr.open_dataset(f'{data_url}rcp85/{target_var}_CORDEX-EUR-11_rcp85_mon_200601-210012_v02.nc')
data_fut_ref = data_fut.sel(time=slice('2081','2100'))
data_fut_ref = data_fut_ref.sel(lat=lat_iberia, lon=lon_iberia)
data_fut_ref_r1 = data_fut_ref.where(data_hist_ref.gcm_variant == 'r1i1p1', drop=True)
data_dx_delta = {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list}
for rcm_name in rcm_list:
    data_fut_rcm = data_fut_ref_r1.where(data_hist_ref_r1.rcm_model==rcm_name, drop=True)
    for gcm_name in gcm_list:
        data_fut_selected = data_fut_rcm.where(data_fut_rcm.gcm_model==gcm_name, drop=True)
        if data_fut_selected.sizes.get("member", 0) == 0:
            continue
        data_fut_selected_mean = data_fut_selected.mean(dim='time').squeeze('member')
        data_dx_delta[rcm_name][gcm_name] = data_fut_selected_mean - data_dx_climatology[rcm_name][gcm_name]

utils.multi_map(data=data_dx_delta, x_map=rcm_dict, y_map=gcm_dict, vlimits=(-10, 10), var='rx1day',
        color='BrBG', cbar_limits=(0, 10, 10), title='Delta rx1day - mm (Average rx1day over 20 years, 1986-2005 as reference, and 2081-2100 as target.)',
        fig_path=FIGURES_PATH, fig_name='delta_rx1day.png')

# DELTA RX1DAY Relative
for rcm_name in rcm_list:
    for gcm_name in gcm_list:
        fut = data_dx_delta[rcm_name][gcm_name]
        ref = data_dx_climatology[rcm_name][gcm_name]
        if fut is not None and ref is not None:
            data_dx_delta[rcm_name][gcm_name] = (fut/ref)*100
        else:
            continue

utils.multi_map(data=data_dx_delta, x_map=rcm_dict, y_map=gcm_dict, vlimits=(-80, 80), var='rx1day',
        color='BrBG', cbar_limits=(0, 10, 10), title='Relative delta rx1day - % (Average rx1day over 20 years, 1986-2005 as reference, and 2081-2100 as target. (Target-Reference/Reference)x100)',
        fig_path=FIGURES_PATH, fig_name='delta_rx1day_relative.png')

# CLIMATOLOGY PRHMAX

lat_min, lat_max = lat_iberia.start, lat_iberia.stop
lon_min, lon_max = lon_iberia.start, lon_iberia.stop

target_years = ['1986', '1991', '1996', '2001']
data_prh_climatology = {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list}

for rcm_name in rcm_list:
    for gcm_name in gcm_list:
        files = [
            f for f in utils.raw_filepath 
            if gcm_name in f and rcm_name in f
        ]
        files_years = [
            f for f in files 
            if any(year in f for year in target_years)
        ]

        if not files_years:
            print(f"⚠️ No files for {gcm_name} + {rcm_name}")
            continue

        data_prh = xr.open_mfdataset(
            files_years, 
            combine="by_coords"
        )
        if 'latitude' in data_prh.coords:
            data_prh_renamed = data_prh.rename({'longitude': 'newlon', 'latitude': 'newlat'})
            lat=data_prh_renamed['newlat'].compute()
            lon=data_prh_renamed['newlon'].compute()
        else:
            data_prh_renamed = data_prh
            lat=data_prh_renamed.lat.compute()
            lon=data_prh_renamed.lon.compute()
        

        data_prh_selected = data_prh.where(
            (lat >= lat_min) & (lat <= lat_max) &
            (lon >= lon_min) & (lon <= lon_max),
            drop=True
        )
        if any(data_prh_selected.sizes.get(dim, 0) == 0 for dim in data_prh_selected.dims):
            print("Alguna de las dimensiones está vacía o no existe")
            continue

        data_prh_max = data_prh_selected[['prhmax']].resample(time="1MS").max()
        data_prh_max_mean = data_prh_max.mean(dim='time')*86400

        data_prh_climatology[rcm_name][gcm_name] = data_prh_max_mean

utils.multi_map(data=data_prh_climatology, x_map=rcm_dict, y_map=gcm_dict, vlimits=(0, 250), var='prhmax',
        color='Blues', cbar_limits=(0, 10, 10), title='Climatology prhmax (montly max prhmax) - mm (Average rx1day over 20 years, 1986-2005)',
        fig_path=FIGURES_PATH, fig_name='climatology_prhmax.png')

# DELTA PRHMAX
data_prh_delta = {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list}
target_years = ['2081', '2086', '2091', '2096']

for rcm_name in rcm_list:
    for gcm_name in gcm_list:
        if data_prh_climatology[rcm_name][gcm_name] is None:
            print("Salimos!")
            continue
        files = [
            f for f in utils.raw_filepath 
            if gcm_name in f and rcm_name in f
        ]
        files_years = [
            f for f in files 
            if any(year in f for year in target_years)
        ]
        files_selected = [
            f for f in files_years
            if 'r1i1p1' in f
        ]
        if not files_selected:
            print(f"⚠️ No files for {gcm_name} + {rcm_name}")
            continue
        data_prh = xr.open_mfdataset(
            files_selected, 
            combine="by_coords"
        )
        if 'latitude' in data_prh.coords:
            data_prh_renamed = data_prh.rename({'longitude': 'newlon', 'latitude': 'newlat'})
            lat=data_prh_renamed['newlat'].compute()
            lon=data_prh_renamed['newlon'].compute()
        else:
            data_prh_renamed = data_prh
            lat=data_prh_renamed.lat.compute()
            lon=data_prh_renamed.lon.compute()
        
        data_prh_selected = data_prh.where(
            (lat >= lat_min) & (lat <= lat_max) &
            (lon >= lon_min) & (lon <= lon_max),
            drop=True
        )
        if any(data_prh_selected.sizes.get(dim, 0) == 0 for dim in data_prh_selected.dims):
            print("Alguna de las dimensiones está vacía o no existe")
            continue

        data_prh_max = data_prh_selected[['prhmax']].resample(time="1MS").max()
        data_prh_max_mean = data_prh_max.mean(dim='time')*86400

        data_prh_delta[rcm_name][gcm_name] = data_prh_max_mean - data_prh_climatology[rcm_name][gcm_name]

utils.multi_map(data=data_prh_delta, x_map=rcm_dict, y_map=gcm_dict, vlimits=(-50, 50), var='prhmax',
        color='BrBG', cbar_limits=(0, 10, 10), title='Delta prhmax (montly max prhmax) - mm (Average rx1day over 20 years, 1986-2005 as reference, and 2081-2100 as target.)',
        fig_path=FIGURES_PATH, fig_name='delta_prhmax.png')

#DELTA PRHMAX Relative
for rcm_name in rcm_list:
    for gcm_name in gcm_list:
        fut = data_prh_delta[rcm_name][gcm_name]
        ref = data_prh_climatology[rcm_name][gcm_name]
        if fut is not None and ref is not None:
            data_prh_delta[rcm_name][gcm_name] = (fut/ref)*100
        else:
            continue

utils.multi_map(data=data_prh_delta, x_map=rcm_dict, y_map=gcm_dict, vlimits=(-80, 80), var='prhmax',
        color='BrBG', cbar_limits=(0, 10, 10), title='Relative delta prhmax (montly max prhmax) - % (Average rx1day over 20 years, 1986-2005 as reference, and 2081-2100 as target.(Target-Reference/Reference)x100))',
        fig_path=FIGURES_PATH, fig_name='delta_prhmax_relative.png')
