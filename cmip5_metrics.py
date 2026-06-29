import xarray as xr
import numpy as np
from pathlib import Path
import sys
import time
import utils_explore as utils
import raw_paths
import os
import glob


DATA_PATH_METRICS = '/nfs/home/gmeteo/reyess/rcm_exploration/data/metrics/'
data_day_hist_url = '/lustre/gmeteo/WORK/DATA/C3S-CDS/C3S-CICA-Atlas/v2/CORDEX-EUR-11/historical/rx1day_CORDEX-EUR-11_historical_mon_197001-200512_v02.nc'
data_day_fut_url = '/lustre/gmeteo/WORK/DATA/C3S-CDS/C3S-CICA-Atlas/v2/CORDEX-EUR-11/rcp85/rx1day_CORDEX-EUR-11_rcp85_mon_200601-210012_v02.nc'

# rcm_name = 'HadREM3-GA7-05'
# gcm_name = 'CNRM-CM5'
# metric = "yr-mean"
# target_var = 'ensemble'
# area = 'Iberia'

target_var = sys.argv[1] # rx1day or prhmax
metric = sys.argv[2] # mon-mean, yr-mean, max-mean
area = sys.argv[3] #Valencia or Iberia

print(f"Doing metric: {metric} - area: {area} - target var: {target_var}")


lat_target, lon_target = utils.get_slice_coords(area)

# Create general structure
# rcm_list = ['ALADIN63', 'CCLM4-8-17', 'COSMO-crCLIM-v1-1', 'HIRHAM5',
#             'HadREM3-GA7-05', 'RACMO22E', 'RCA4', 'REMO2009',
#             'REMO2015', 'RegCM4-6', 'WTF361H', 'WRF381P']
# gcm_list = ['CNRM-CM5', 'CanESM2', 'EC-EARTH', 'HadGEM2-ES',
#             'IPSL-CM5A-MR', 'MIROC5', 'MPI-ESM-LR', 'NorESM1-M']
rcm_list = ['CCLM4-8-17', 'COSMO-crCLIM-v1-1', 'HIRHAM5',
            'HadREM3-GA7-05', 'RACMO22E']
gcm_list = ['CNRM-CM5', 'EC-EARTH', 'HadGEM2-ES', 'MPI-ESM-LR', 'NorESM1-M']
seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']

rcm_dict = {rcm_name:rcm_name for rcm_name in rcm_list}
gcm_dict = {gcm_name:gcm_name for gcm_name in gcm_list}

prhmax_paths = raw_paths.raw_filepath_prhmax
gcm_gwl3_years = utils.gcm_gwl3_years
time_start = time.time()

if target_var == 'rx1day':
    # CLIMATOLOGY DX1DAY
    data_hist = xr.open_dataset(data_day_hist_url)
    data_hist_ref = data_hist.sel(time=slice('1986','2005'))
    data_hist_ref = data_hist_ref.sel(lat=lat_target, lon=lon_target)#IBERIA
    for rcm_name in rcm_list:
        data_hist_rcm = data_hist_ref.where(data_hist_ref.rcm_model==rcm_name, drop=True)
        for gcm_name in gcm_list:
            data_hist_selected = data_hist_rcm.where(data_hist_rcm.gcm_model==gcm_name, drop=True)
            
            if data_hist_selected.sizes.get("member", 0) == 0:
                continue
            first_variant = data_hist_selected.gcm_variant.values[0]

            data_hist_selected = data_hist_selected.where(
                data_hist_selected.gcm_variant == first_variant,
                drop=True
            )
            data_selected = data_hist_selected.squeeze('member')
            season_result = []
            for season in seasons:
                if season == 'ANN':
                    data_season = data_selected
                else:
                    data_season = data_selected.where(
                        data_selected['time.season'] == season, drop=True
                    )
                if metric == 'mon-mean':
                    data_out = data_season.rx1day.mean(dim='time')
                elif metric == 'yr-mean':
                    data_out = data_season.rx1day.resample(time="YS").max().mean(dim='time')
                elif metric == 'max-mean':
                    data_out = data_season.rx1day.max(dim='time')
                # Añadir coordenada season
                data_out = data_out.expand_dims(season=[season])
                season_result.append(data_out)

            data_final = xr.concat(season_result, dim='season')
            data_final.to_netcdf(f'{DATA_PATH_METRICS}/hist_climatology_{area}_{metric}_{target_var}_{gcm_name}_{rcm_name}_1986-2005.nc')


    # GWL3 DX1DAY
    data_fut = xr.open_dataset(data_day_fut_url)
    data_fut_ref = data_fut.sel(lat=lat_target, lon=lon_target)
    for rcm_name in rcm_list:
        data_fut_rcm = data_fut_ref.where(data_fut_ref.rcm_model==rcm_name, drop=True)
        for gcm_name in gcm_list:
            data_fut_selected = data_fut_rcm.where(data_fut_rcm.gcm_model==gcm_name, drop=True)
            
            data_fut_selected = data_fut_selected.sel(time=slice(gcm_gwl3_years[gcm_name][0], gcm_gwl3_years[gcm_name][1]))
            if data_fut_selected.sizes.get("member", 0) == 0:
                continue
            first_variant = data_fut_selected.gcm_variant.values[0]
            data_fut_selected = data_fut_selected.where(
                data_fut_selected.gcm_variant == first_variant,
                drop=True
            )
            data_selected = data_fut_selected.squeeze('member')
            season_result = []
            for season in seasons:
                if season == 'ANN':
                    data_season = data_selected
                else:
                    data_season = data_selected.where(
                        data_selected['time.season'] == season, drop=True
                    )
                if metric == 'mon-mean':
                    data_out = data_season.rx1day.mean(dim='time')
                elif metric == 'yr-mean':
                    data_out = data_season.rx1day.resample(time="YS").max().mean(dim='time')
                elif metric == 'max-mean':
                    data_out = data_season.rx1day.max(dim='time')

                # Añadir coordenada season
                data_out = data_out.expand_dims(season=[season])
                season_result.append(data_out)

            data_final = xr.concat(season_result, dim='season')

            data_final.to_netcdf(f'{DATA_PATH_METRICS}/gwl3_climatology_{area}_{metric}_{target_var}_{gcm_name}_{rcm_name}_{gcm_gwl3_years[gcm_name][0]}-{gcm_gwl3_years[gcm_name][1]}.nc')


elif target_var == 'prhmax':
# CLIMATOLOGY PRHMAX
    lat_min, lat_max = lat_target.start, lat_target.stop
    lon_min, lon_max = lon_target.start, lon_target.stop

    target_years = ['1986', '1991', '1996', '2001']

    for rcm_name in rcm_list:
        for gcm_name in gcm_list:
            files_selected = utils.filter_prhmax_files(prhmax_paths, gcm_name, rcm_name, target_years)

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


            data_prh_selected = data_prh_renamed.where(
                (lat >= lat_min) & (lat <= lat_max) &
                (lon >= lon_min) & (lon <= lon_max),
                drop=True
            )
#             data_prh_selected2 = data_prh_renamed2.where(
#                 (lat >= lat_min) & (lat <= lat_max) &
#                 (lon >= lon_min) & (lon <= lon_max),
#                 drop=True
#             )
#             # 1. Creamos una máscara booleana 2D limpia (solo rlat y rlon) usando las coords del dataset
# mask_2d = (
#     (data_prh_renamed['lat'] >= lat_min) & (data_prh_renamed['lat'] <= lat_max) &
#     (data_prh_renamed['lon'] >= lon_min) & (data_prh_renamed['lon'] <= lon_max)
# )

# # 2. Encontramos qué índices de las dimensiones rlat y rlon contienen a Valencia
# valid_rlat = mask_2d.any(dim='rlon')
# valid_rlon = mask_2d.any(dim='rlat')

# # 3. Hacemos un recorte físico estricto de la matriz pasándole los vectores indexados.
# # Esto reduce el dataset de 46GB a solo unos pocos Megabytes en un instante.
# data_prh_selected3 = data_prh_renamed.sel(
#     rlat=valid_rlat, 
#     rlon=valid_rlon
# ).where(mask_2d) # Aplica la máscara final para limpiar los bordes fuera de la caja

# # 4. Ahora cuenta los NaNs de una sola tajada temporal para verificar que no esté vacío
# # Seleccionamos solo el primer día (.isel(time=0)) para no cargar los 46GB a la RAM
# nans_un_dia = data_prh_selected2['prhmax'].load().isnull().sum().item()
# nans_selected2 = data_prh_selected2['prhmax'].load().isnull().sum().item()

# print(f"NaNs en el primer mapa recortado: {nans_un_dia}")
            if any(data_prh_selected.sizes.get(dim, 0) == 0 for dim in data_prh_selected.dims):
                print("Alguna de las dimensiones está vacía o no existe")
                continue

            season_result = []
            for season in seasons:
                if season == 'ANN':
                    data_season = data_prh_selected
                else:
                    data_season = data_prh_selected.where(
                        data_prh_selected['time.season'] == season, drop=True
                    )
                if metric == 'mon-mean':
                    data_out = data_season[['prhmax']].resample(time="1MS").max()
                    data_out = data_out.mean(dim='time')*86400
                elif metric == 'yr-mean':
                    data_out = data_season[['prhmax']].resample(time="YS").max()
                    data_out = data_out.mean(dim='time')*86400
                elif metric == 'max-mean':
                    data_out = data_season[['prhmax']].max(dim='time')*86400

                # Añadir coordenada season
                data_out = data_out.expand_dims(season=[season])
                season_result.append(data_out)

            data_final = xr.concat(season_result, dim='season') 
            data_final.to_netcdf(f'{DATA_PATH_METRICS}/hist_climatology_{area}_{metric}_{target_var}_{gcm_name}_{rcm_name}_1986-2005.nc')

    # GWL3 PRHMAX
    for rcm_name in rcm_list:
        for gcm_name in gcm_list:
            files_selected = utils.filter_prhmax_files(prhmax_paths, gcm_name, rcm_name)
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

            season_result = []
            for season in seasons:
                if season == 'ANN':
                    data_season = data_prh_selected
                else:
                    data_season = data_prh_selected.where(
                        data_prh_selected['time.season'] == season, drop=True
                    )
                if metric == 'mon-mean':
                    data_out = data_season[['prhmax']].resample(time="1MS").max()
                    data_out = data_out.mean(dim='time')*86400
                elif metric == 'yr-mean':
                    data_out = data_season[['prhmax']].resample(time="YS").max()
                    data_out = data_out.mean(dim='time')*86400
                elif metric == 'max-mean':
                    data_out = data_season[['prhmax']].max(dim='time')*86400

                # Añadir coordenada season
                data_out = data_out.expand_dims(season=[season])
                season_result.append(data_out)

            data_final = xr.concat(season_result, dim='season') 

            data_final.to_netcdf(f'{DATA_PATH_METRICS}/gwl3_climatology_{area}_{metric}_{target_var}_{gcm_name}_{rcm_name}_{gcm_gwl3_years[gcm_name][0]}-{gcm_gwl3_years[gcm_name][1]}.nc')
elif target_var == 'ensemble':
    # Here are the mean, p20 and p80 for the whole ensemble of rx1day and prhmax
    rx1day_list = []
    prhmax_list = []
    prh_ref = xr.open_dataset(f'{DATA_PATH_METRICS}/hist_climatology_{area}_{metric}_prhmax_CNRM-CM5_CCLM4-8-17_1986-2005.nc')
    for rcm_name in rcm_list:
        for gcm_name in gcm_list:
            # Load metrics
            path_hist_rx1day = f'{DATA_PATH_METRICS}/hist_climatology_{area}_{metric}_rx1day_{gcm_name}_{rcm_name}_1986-2005.nc'
            path_fut_rx1day = f'{DATA_PATH_METRICS}/gwl3_climatology_{area}_{metric}_rx1day_{gcm_name}_{rcm_name}_{gcm_gwl3_years[gcm_name][0]}-{gcm_gwl3_years[gcm_name][1]}.nc'
            if not os.path.exists(path_hist_rx1day) or not os.path.exists(path_fut_rx1day):
                print(f"⚠️ Missing files for {gcm_name} + {rcm_name} in rx1day")
                continue
            print(f"Processing RX1DAY {gcm_name} + {rcm_name}")
            metric_hist_rx1day = xr.open_dataset(path_hist_rx1day)
            metric_fut_rx1day = xr.open_dataset(path_fut_rx1day)
            relative_rx1day = (metric_fut_rx1day - metric_hist_rx1day) / metric_hist_rx1day * 100
            rx1day_list.append(relative_rx1day)
            

            path_hist_prhmax = f'{DATA_PATH_METRICS}/hist_climatology_{area}_{metric}_prhmax_{gcm_name}_{rcm_name}_1986-2005.nc'
            path_fut_prhmax = f'{DATA_PATH_METRICS}/gwl3_climatology_{area}_{metric}_prhmax_{gcm_name}_{rcm_name}_{gcm_gwl3_years[gcm_name][0]}-{gcm_gwl3_years[gcm_name][1]}.nc'
            if not os.path.exists(path_hist_prhmax) or not os.path.exists(path_fut_prhmax):
                print(f"⚠️ Missing files for {gcm_name} + {rcm_name} in prhmax")
                continue
            print(f"Processing PRHMAX {gcm_name} + {rcm_name}")
            metric_hist_prhmax = xr.open_dataset(path_hist_prhmax)
            coords = {"lat": "latitude", "lon": "longitude"} if rcm_name == 'HadREM3-GA7-05' else {"lat": "lat", "lon": "lon"}
            metric_hist_prhmax = utils.smart_regrid(metric_hist_prhmax, prh_ref, var="prhmax", space_coords=coords)
            metric_hist_prhmax = utils.fix_latlon(metric_hist_prhmax)
            metric_fut_prhmax = xr.open_dataset(path_fut_prhmax)
            metric_fut_prhmax = utils.smart_regrid(metric_fut_prhmax, prh_ref, var="prhmax", space_coords=coords)
            metric_fut_prhmax = utils.fix_latlon(metric_fut_prhmax)

            relative_prhmax = (metric_fut_prhmax - metric_hist_prhmax) / metric_hist_prhmax * 100
            relative_prhmax = relative_prhmax.assign_coords({
                'rlat': np.round(relative_prhmax.rlat.values, 4),
                'rlon': np.round(relative_prhmax.rlon.values, 4)
            })
            prhmax_list.append(relative_prhmax)
    rx1day_concat = xr.concat(rx1day_list, dim='model')
    rx1day_mean = rx1day_concat.mean(dim='model')
    rx1day_p20 = rx1day_concat.quantile(0.2, dim='model')
    rx1day_p80 = rx1day_concat.quantile(0.8, dim='model')
    rx1day_mean.to_netcdf(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_mean_{area}_{metric}_rx1day.nc')
    rx1day_p20.to_netcdf(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_p20_{area}_{metric}_rx1day.nc')
    rx1day_p80.to_netcdf(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_p80_{area}_{metric}_rx1day.nc')

    prhmax_concat = xr.concat(prhmax_list, dim='model')
    prhmax_mean = prhmax_concat.mean(dim='model')
    prhmax_p20 = prhmax_concat.quantile(0.2, dim='model')
    prhmax_p80 = prhmax_concat.quantile(0.8, dim='model')
    prhmax_mean.to_netcdf(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_mean_{area}_{metric}_prhmax.nc')
    prhmax_p20.to_netcdf(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_p20_{area}_{metric}_prhmax.nc')
    prhmax_p80.to_netcdf(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_p80_{area}_{metric}_prhmax.nc')


time_end = time.time()
time_elapsed = time_end - time_start
print(f"Time elapsed: {time_elapsed/60:.2f} minutes for {target_var} and {metric}")


# import os
# import matplotlib.pyplot as plt

# # Crear la figura
# fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# # 1. Graficar RX1DAY
# metric_hist_rx1day['rx1day'].sel(season='ANN').plot(
#     ax=axes[0], 
#     cmap='Blues', 
#     cbar_kwargs={'label': 'mm'}
# )
# axes[0].set_title('Estructura de RX1DAY\n(Rejilla Estándar Lat/Lon)', fontweight='bold')

# # 2. Graficar PRHMAX
# metric_hist_prhmax['prhmax'].sel(season='ANN').plot(
#     ax=axes[1], 
#     x='rlon', y='rlat', 
#     cmap='Reds', 
#     cbar_kwargs={'label': 'Valor'}
# )
# axes[1].set_title('Estructura de PRHMAX\n(Rejilla Rotada rlat/rlon)', fontweight='bold')

# plt.tight_layout()

# # --- GUARDAR EN TU RUTA ---
# output_dir = '/nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration/example_figures'
# os.makedirs(output_dir, exist_ok=True)  # Crea la carpeta si no existe

# plt.savefig(f'{output_dir}/comparativa_estructuras.png', bbox_inches='tight', dpi=150)
# plt.close()