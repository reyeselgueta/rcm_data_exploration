import time
import xarray as xr
import numpy as np
from pathlib import Path
import sys
import utils_explore as utils

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

target_coords = sys.argv[1] 
plot_num = sys.argv[2] # '1' or '2'
flux = sys.argv[3] # 1 or 0
custom = int(sys.argv[4]) # 1, 2, or 3
# statistic = sys.argv[5] # 'mean' 'or 'max' or 'date_max'

# target_coords = 'valencia'
# plot_num = '6'
# flux = False
# custom = 1
print(f"Target coords: {target_coords} Plot num: {plot_num} Flux: {flux} Custom: {custom}")
#

# EURR-3 solo tiene historical y de ICTP

if target_coords == 'iberia':
    lat_center = 39.000
    lon_center = -2.875
    folder_name = 'iberia-3'
    grid_multiplier = 4
elif target_coords == 'valencia':
    lat_center = 39.467
    lon_center = -0.375
    folder_name = 'val-3'
    grid_multiplier = 1
elif target_coords == 'europe':
    lat_center = 39.467
    lon_center = -0.375
    folder_name = 'eur-3'
    grid_multiplier = 12
else:
    raise ValueError("target_coords must be 'iberia', 'valencia', or 'europe'")


data_dict = utils.data_rcm_paths

seasons = ['Annual', 'DJF', 'MAM', 'JJA', 'SON']
frecuencies = ['day', '1hr']
objective_rcms = ['CNRM-MF', 'BCCR-UCAN', 'BCCR-UCAN_eur12']

time1= time.time()
if plot_num == '1':
    for frecuency in frecuencies: 
        for rcm_name, data_rcm in data_dict.items():
            # General parameters
            gcm_name = gcm_name
            rcm_dict = {rcm_name:rcm_name}
            gcm_dict = {gcm_name:gcm_name}
            time2 = time.time()
            time_multiplier = 3600 if frecuency == '1hr' else 86400
            value_limits = {
                'climatology': (0,720) if frecuency == 'day' else (0, 180),
                'delta': (-40,40) if frecuency == 'day' else (-20, 20),
                'relative': (-100, 100)
                } 

            #FOR CLIMATOLOGY
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/climatology_{frecuency}_{rcm_name}_{target_coords}_1995-2015.nc").compute().pr
            data_hist = {season: {rcm_name : {gcm_name: gcm_name}} for season in seasons}
            for season in seasons:
                data_hist[season][rcm_name][gcm_name] = data.sel(season=season)
                utils.multi_map(data=data_hist[season], x_map=rcm_dict, y_map=gcm_dict, vlimits=value_limits['climatology'],# var='prhmax',
                        custom_color=utils.custom_color_2, cbar_limits=(0, 20, 20), title=f'Climatology mm/{frecuency}-{season}',
                        fig_path=FIGURES_PATH, fig_name=f'Climatology_{frecuency}_{season}_{rcm_name}_{target_coords}.png')

            time3 = time.time()
            print(f"Datos guardados hist {frecuency} {rcm_name}: {(time3-time2)/60:.2f} minutes")

            #FOR SSP
            if 'ssp370' not in data_rcm.keys():
                continue
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/future_{frecuency}_{rcm_name}_{target_coords}_CWL3.nc").compute().pr
            data_delta = {season: {rcm_name : {gcm_name: gcm_name}} for season in seasons}
            data_future = {season: {rcm_name : {gcm_name: gcm_name}} for season in seasons}
            for season in seasons:
                data_delta[season][rcm_name][gcm_name] = data.sel(season=season) - data_hist[season][rcm_name][gcm_name]
                data_future[season][rcm_name][gcm_name] = data.sel(season=season)
                utils.multi_map(data=data_delta[season], x_map=rcm_dict, y_map=gcm_dict, vlimits=value_limits['delta'],# var='prhmax',
                        color='BrBG', cbar_limits=(0, 20, 20), title=f'Delta mm/{frecuency}-{season}',
                        fig_path=FIGURES_PATH, fig_name=f'Delta_{frecuency}_{season}_{rcm_name}_{target_coords}.png')
                utils.multi_map(data=data_future[season], x_map=rcm_dict, y_map=gcm_dict, vlimits=value_limits['delta'],# var='prhmax',
                        custom_color=utils.custom_color_2, cbar_limits=(0, 20, 20), title=f'Delta mm/{frecuency}-{season}',
                        fig_path=FIGURES_PATH, fig_name=f'Future_{frecuency}_{season}_{rcm_name}_{target_coords}.png')

            time4= time.time()
            print(f"Datos guardados future {frecuency} {rcm_name}: {(time4-time3)/60:.2f} minutes")

            # RELATIVE
            data_relative = {season: {rcm_name : {gcm_name: gcm_name}} for season in seasons}
            for season in seasons:
                data_relative[season][rcm_name][gcm_name] = (data_delta[season][rcm_name][gcm_name]/data_hist[season][rcm_name][gcm_name])*100
                utils.multi_map(data=data_relative[season], x_map=rcm_dict, y_map=gcm_dict, vlimits=value_limits['relative'],# var='prhmax',
                        color='BrBG', cbar_limits=(0, 20, 20), title=f'Relative mm/{frecuency}-{season}',
                        fig_path=FIGURES_PATH, fig_name=f'Relative_{frecuency}_{season}_{rcm_name}_{target_coords}.png')



    time5 = time.time()
    print(f"Plot terminados en: {(time5-time1)/60:.2f} minutes")


# for rcm_name in objective_rcms:
#     for frecuency in frecuencies:
#         print(f"{rcm_name} {frecuency} {target_coords}")
#         max_value = data_ssp370['Annual'][rcm_name][frecuency].max().compute() if data_ssp370['Annual'][rcm_name][frecuency] is not None else None
#         print(f"Max future {max_value.values if max_value is not None else None} Over 50:{max_value.values > 50 if max_value is not None else None}")
#         max_value = data_hist['Annual'][rcm_name][frecuency].max().compute() if data_hist['Annual'][rcm_name][frecuency] is not None else None
#         print(f"Max hist {max_value.values if max_value is not None else None} Over 50:{max_value.values > 50 if max_value is not None else None}")

elif plot_num == '2': # Mean Max 2 RCMs
    
    objective_rcms = ['BCCR-UCAN', 'BCCR-UCAN_eur12']
    rcm_dict = {'BCCR-UCAN': '3 Kilometers', 'BCCR-UCAN_eur12': '12 Kilometers'}
    mask_dict = {'BCCR-UCAN': None, 'BCCR-UCAN_eur12': None}
    frecuency_dict = {'day': 'Daily', '1hr': 'Hourly'}
    data_hist = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_ssp370 = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_delta = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_relative = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    for frecuency in frecuencies: 
        if flux == 1:
            divider = 3600 if frecuency == '1hr' else 86400
        else:
            divider = 1
        for rcm_name in objective_rcms:
            gcm_name = data_dict[rcm_name]['gcm']
            
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/climatology_{frecuency}_{rcm_name}_{target_coords}_1995-2015_mask.nc").compute().pr
            for season in seasons:
                data_hist[season][rcm_name][frecuency] = data.sel(season=season)/divider

            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/ssp370_{frecuency}_{rcm_name}_{target_coords}_CWL3_mask.nc").compute().pr
            for season in seasons:
                data_ssp370[season][rcm_name][frecuency] = data.sel(season=season)/divider
            
            for season in seasons:
                data = data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]
                data_delta[season][rcm_name][frecuency] = data

            for season in seasons:
                data = ((data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]) / data_hist[season][rcm_name][frecuency]) * 100
                data_relative[season][rcm_name][frecuency] = data
        
            del data

            mask_data = xr.open_mfdataset(f"{data_dict[rcm_name]['ssp370']['fx']}*.nc")
            lat=mask_data.lat.values
            lon=mask_data.lon.values
            abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
            iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

            # Slice area around the city
            grid_number = 90 if rcm_name != 'BCCR-UCAN_eur12' else 30
            grid_number = grid_number * grid_multiplier
            half_box = grid_number // 2
            y_slice = slice(max(0, iy - half_box), iy + half_box)
            x_slice = slice(max(0, ix - half_box), ix + half_box)

            mask = mask_data['sftlf'].isel(y=y_slice, x=x_slice)
            mask_dict[rcm_name] = mask
    

    if flux == 1:
        value_limits = {
            'climatology': [(0,0.015), (0, 0.015)],
            'delta': [(-0.008,0.008), (-0.008, 0.0080)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        custom_ticks = [utils.custom_ticks_flux, utils.custom_ticks_flux]
    else:
        value_limits = {
            'climatology': [(0,200), (0, 50)],
            'delta': [(-40,40), (-20, 20)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        if custom == 1:
            custom_colorbar = [utils.custom_color_1, utils.custom_color_1]
            color_extra_name = '_customcolor1'
            custom_ticks = None
        elif custom == 2:
            custom_colorbar = utils.custom_color_2
            color_extra_name = '_customcolor2'
            custom_ticks =  [utils.custom_ticks_2_day, utils.custom_ticks_2_1hr]
        elif custom == 3:
            custom_colorbar = [utils.custom_color_3, utils.custom_color_3]
            color_extra_name = '_customcolor3'
            custom_ticks = [utils.custom_ticks_3, utils.custom_ticks_3]

    for season in seasons:
        utils.multi_map(data=data_hist[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'Historical-{season} BCCR-UCAN',
                fig_path=FIGURES_PATH, fig_name=f'Climatology_{season}_BCCR-UCAN_vs_BCCR-UCAN_eur12_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_ssp370[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'SSP370-{season}  BCCR-UCAN',
                fig_path=FIGURES_PATH, fig_name=f'Ssp370_{season}_BCCR-UCAN_vs_BCCR-UCAN_eur12_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_delta[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['delta'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['delta'],
                title=f'Delta-{season}  BCCR-UCAN',
                fig_path=FIGURES_PATH, fig_name=f'Delta_{season}_BCCR-UCAN_vs_BCCR-UCAN_eur12_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_relative[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['relative'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['relative'],
                title=f'Relative-{season} BCCR-UCAN',
                fig_path=FIGURES_PATH, fig_name=f'Relative_{season}_BCCR-UCAN_vs_BCCR-UCAN_eur12_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)
        
    print("Plots terminados")

elif plot_num == '3': # Mean Max CNRM
    
    objective_rcms = ['CNRM-MF']
    rcm_dict = {'CNRM-MF': '3 Kilometers'}
    mask_dict = {'CNRM-MF': None}
    frecuency_dict = {'day': 'Daily', '1hr': 'Hourly'}
    data_hist = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_ssp370 = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_delta = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_relative = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    for frecuency in frecuencies: 
        if flux == 1:
            divider = 3600 if frecuency == '1hr' else 86400
        else:
            divider = 1
        for rcm_name in objective_rcms:
            gcm_name = data_dict[rcm_name]['gcm']
            
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/climatology_{frecuency}_{rcm_name}_{target_coords}_1995-2015_mask.nc").compute().pr
            for season in seasons:
                data_hist[season][rcm_name][frecuency] = data.sel(season=season)/divider

            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/ssp370_{frecuency}_{rcm_name}_{target_coords}_CWL3_mask.nc").compute().pr
            for season in seasons:
                data_ssp370[season][rcm_name][frecuency] = data.sel(season=season)/divider
            
            for season in seasons:
                data = data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]
                data_delta[season][rcm_name][frecuency] = data

            for season in seasons:
                data = ((data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]) / data_hist[season][rcm_name][frecuency]) * 100
                data_relative[season][rcm_name][frecuency] = data
        
            del data

            mask_data = xr.open_mfdataset(f"{data_dict[rcm_name]['ssp370']['fx']}*.nc")
            lat=mask_data.lat.values
            lon=mask_data.lon.values
            abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
            iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

            # Slice area around the city
            grid_number = 90 if rcm_name != 'BCCR-UCAN_eur12' else 30
            grid_number = grid_number * grid_multiplier
            half_box = grid_number // 2
            y_slice = slice(max(0, iy - half_box), iy + half_box)
            x_slice = slice(max(0, ix - half_box), ix + half_box)

            mask = mask_data['sftlf'].isel(y=y_slice, x=x_slice)
            mask_dict[rcm_name] = mask
    
    if flux == 1:
        value_limits = {
            'climatology': [(0,0.015), (0, 0.015)],
            'delta': [(-0.008,0.008), (-0.008, 0.0080)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        custom_ticks = [utils.custom_ticks_flux, utils.custom_ticks_flux]
    else:
        value_limits = {
            'climatology': [(0,200), (0, 50)],
            'delta': [(-40,40), (-20, 20)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        if custom == 1:
            custom_colorbar = [utils.custom_color_1, utils.custom_color_1]
            color_extra_name = '_customcolor1'
            custom_ticks = None
        elif custom == 2:
            custom_colorbar = [utils.custom_color_2, utils.custom_color_2]
            color_extra_name = '_customcolor2'
            custom_ticks =  [utils.custom_ticks_2_day, utils.custom_ticks_2_1hr]
        elif custom == 3:
            custom_colorbar = [utils.custom_color_3, utils.custom_color_3]
            color_extra_name = '_customcolor3'
            custom_ticks = [utils.custom_ticks_3, utils.custom_ticks_3]

    for season in seasons:
        utils.multi_map(data=data_hist[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'Historical-{season} CNRM-MF',
                fig_path=FIGURES_PATH, fig_name=f'Climatology_{season}_CNRM-MF_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_ssp370[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'SSP370-{season}  CNRM-MF',
                fig_path=FIGURES_PATH, fig_name=f'Ssp370_{season}_CNRM-MF_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_delta[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['delta'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['delta'],
                title=f'Delta-{season}  CNRM-MF',
                fig_path=FIGURES_PATH, fig_name=f'Delta_{season}_CNRM-MF_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_relative[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['relative'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['relative'],
                title=f'Relative-{season} CNRM-MF',
                fig_path=FIGURES_PATH, fig_name=f'Relative_{season}_CNRM-MF_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)
        
        
    print("Plots terminados")

elif plot_num == '4': # Mean Max 3 RCMs
    
    objective_rcms = ['CNRM-MF', 'BCCR-UCAN', 'BCCR-UCAN_eur12']
    rcm_dict = {'CNRM-MF': '3 Kilometers', 'BCCR-UCAN': '3 Kilometers', 'BCCR-UCAN_eur12': '12 Kilometers'}
    mask_dict = {'CNRM-MF': None, 'BCCR-UCAN': None, 'BCCR-UCAN_eur12': None}
    frecuency_dict = {'day': 'Daily', '1hr': 'Hourly'}
    data_hist = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_ssp370 = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_delta = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_relative = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    for frecuency in frecuencies: 
        if flux == 1:
            divider = 3600 if frecuency == '1hr' else 86400
        else:
            divider = 1
        for rcm_name in objective_rcms:
            gcm_name = data_dict[rcm_name]['gcm']
            
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/climatology_{frecuency}_{rcm_name}_{target_coords}_mean_1995-2015.nc").compute().pr
            for season in seasons:
                data_hist[season][rcm_name][frecuency] = data.sel(season=season)/divider

            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/ssp370_{frecuency}_{rcm_name}_{target_coords}_mean_CWL3.nc").compute().pr
            for season in seasons:
                data_ssp370[season][rcm_name][frecuency] = data.sel(season=season)/divider
            
            for season in seasons:
                data = data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]
                data_delta[season][rcm_name][frecuency] = data

            for season in seasons:
                data = ((data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]) / data_hist[season][rcm_name][frecuency]) * 100
                data_relative[season][rcm_name][frecuency] = data
        
            del data

            mask_data = xr.open_mfdataset(f"{data_dict[rcm_name]['ssp370']['fx']}*.nc")
            lat=mask_data.lat.values
            lon=mask_data.lon.values
            abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
            iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

            # Slice area around the city
            grid_number = 90 if rcm_name != 'BCCR-UCAN_eur12' else 30
            grid_number = grid_number * grid_multiplier
            half_box = grid_number // 2
            y_slice = slice(max(0, iy - half_box), iy + half_box)
            x_slice = slice(max(0, ix - half_box), ix + half_box)

            mask = mask_data['sftlf'].isel(y=y_slice, x=x_slice)
            mask_dict[rcm_name] = mask
    
    if flux == 1:
        value_limits = {
            'climatology': [(0,0.015), (0, 0.015)],
            'delta': [(-0.008,0.008), (-0.008, 0.0080)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        custom_ticks = [utils.custom_ticks_flux, utils.custom_ticks_flux]
    else:
        value_limits = {
            'climatology': [(0,120), (0, 40)],
            'delta': [(-40,40), (-20, 20)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        if custom == 1:
            custom_colorbar = [utils.custom_color_1, utils.custom_color_1]
            color_extra_name = '_customcolor1'
            custom_ticks = None
        elif custom == 2:
            custom_colorbar = [utils.custom_color_2, utils.custom_color_2]
            color_extra_name = '_customcolor2'
            custom_ticks =  [utils.custom_ticks_2_day, utils.custom_ticks_2_1hr]
        elif custom == 3:
            custom_colorbar = [utils.custom_color_3, utils.custom_color_3]
            color_extra_name = '_customcolor3'
            custom_ticks = [utils.custom_ticks_3, utils.custom_ticks_3]

    for season in seasons:
        utils.multi_map(data=data_hist[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'Historical-{season} CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Climatology_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_ssp370[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'SSP370-{season}  CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Ssp370_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_delta[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['delta'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['delta'],
                title=f'Delta-{season}  CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Delta_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_relative[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['relative'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['relative'],
                title=f'Relative-{season} CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Relative_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)
        
    print("Plots terminados")


elif plot_num == '5': #MAX MAX
    
    objective_rcms = ['CNRM-MF', 'BCCR-UCAN', 'BCCR-UCAN_eur12']
    rcm_dict = {'CNRM-MF': '3 Kilometers', 'BCCR-UCAN': '3 Kilometers', 'BCCR-UCAN_eur12': '12 Kilometers'}
    mask_dict = {'CNRM-MF': None, 'BCCR-UCAN': None, 'BCCR-UCAN_eur12': None}
    frecuency_dict = {'day': 'Daily', '1hr': 'Hourly'}
    data_hist = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_ssp370 = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_delta = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_relative = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    for frecuency in frecuencies: 
        if flux == 1:
            divider = 3600 if frecuency == '1hr' else 86400
        else:
            divider = 1
        for rcm_name in objective_rcms:
            gcm_name = data_dict[rcm_name]['gcm']
            
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/climatology_{frecuency}_{rcm_name}_{target_coords}_max_1995-2015.nc").compute().pr
            for season in seasons:
                data_hist[season][rcm_name][frecuency] = data.sel(season=season)/divider

            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/ssp370_{frecuency}_{rcm_name}_{target_coords}_max_CWL3.nc").compute().pr
            for season in seasons:
                data_ssp370[season][rcm_name][frecuency] = data.sel(season=season)/divider
            
            for season in seasons:
                data = data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]
                data_delta[season][rcm_name][frecuency] = data

            for season in seasons:
                data = ((data_ssp370[season][rcm_name][frecuency] - data_hist[season][rcm_name][frecuency]) / data_hist[season][rcm_name][frecuency]) * 100
                data_relative[season][rcm_name][frecuency] = data
        
            del data

            mask_data = xr.open_mfdataset(f"{data_dict[rcm_name]['ssp370']['fx']}*.nc")
            lat=mask_data.lat.values
            lon=mask_data.lon.values
            abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
            iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

            # Slice area around the city
            grid_number = 90 if rcm_name != 'BCCR-UCAN_eur12' else 30
            grid_number = grid_number * grid_multiplier
            half_box = grid_number // 2
            y_slice = slice(max(0, iy - half_box), iy + half_box)
            x_slice = slice(max(0, ix - half_box), ix + half_box)

            mask = mask_data['sftlf'].isel(y=y_slice, x=x_slice)
            mask_dict[rcm_name] = mask
    
    if flux == 1:
        value_limits = {
            'climatology': [(0,0.015), (0, 0.015)],
            'delta': [(-0.008,0.008), (-0.008, 0.0080)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        custom_ticks = [utils.custom_ticks_flux, utils.custom_ticks_flux]
    else:
        value_limits = {
            'climatology': [(0,540), (0, 180)],
            'delta': [(-150,150), (-50, 50)],
            'relative': [(-100, 100), [-100, 100]]
        }

        colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
        if custom == 1:
            custom_colorbar = [utils.custom_color_1, utils.custom_color_1]
            color_extra_name = '_customcolor1'
            custom_ticks = None
        elif custom == 2:
            custom_colorbar = [utils.custom_color_2, utils.custom_color_2]
            color_extra_name = '_customcolor2'
            custom_ticks =  [utils.custom_ticks_2_day, utils.custom_ticks_2_1hr]
        elif custom == 3:
            custom_colorbar = [utils.custom_color_3, utils.custom_color_3]
            color_extra_name = '_customcolor3'
            custom_ticks = [utils.custom_ticks_3, utils.custom_ticks_3]

    for season in seasons:
        utils.multi_map(data=data_hist[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'Historical Max-{season} CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Climatology_Max_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_ssp370[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['climatology'],
                custom_color=custom_colorbar, cbar_limits=colobar_limits['climatology'], custom_cbar_ticks=custom_ticks,
                title=f'SSP370 Max-{season}  CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Ssp370_Max_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_delta[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['delta'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['delta'],
                title=f'Delta Max-{season}  CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Delta_Max_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=data_relative[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits['relative'],
                color=['BrBG', 'BrBG'], cbar_limits=colobar_limits['relative'],
                title=f'Relative Max-{season} CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Relative_Max_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)
        
    print("Plots terminados")


elif plot_num == '6':
    frecuencies = ['day', '1hr']  # Only daily for this plot
    rcm_dict = {rcm_name: rcm_name for rcm_name in objective_rcms}
    mask_dict = {rcm_name: None for rcm_name in objective_rcms}
    frecuency_dict = {'Percentage': 'Percentage'}
    data_hist = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    data_ssp370 = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    percentage_matches_hist = {season: {rcm_name: {'Percentage': None} for rcm_name in objective_rcms} for season in seasons}
    percentage_matches_ssp370 = {season: {rcm_name: {'Percentage': None} for rcm_name in objective_rcms} for season in seasons}
    data_relative = {season: {rcm_name: {frecuency: None for frecuency in frecuencies} for rcm_name in objective_rcms} for season in seasons}
    for frecuency in frecuencies: 
        if flux == 1:
            divider = 3600 if frecuency == '1hr' else 86400
        else:
            divider = 1
        for rcm_name in objective_rcms:
            gcm_name = data_dict[rcm_name]['gcm']
            
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/climatology_{frecuency}_{rcm_name}_{target_coords}_date_max_1995-2015.nc")
            for season in seasons:
                data_hist[season][rcm_name][frecuency] = data.sel(season=season).compute().date_max
            data = xr.open_dataset(f"{DATA_PATH_METRICS}{folder_name}/ssp370_{frecuency}_{rcm_name}_{target_coords}_date_max_CWL3.nc")
            for season in seasons:
                data_ssp370[season][rcm_name][frecuency] = data.sel(season=season).compute().date_max
            
            del data

            mask_data = xr.open_mfdataset(f"{data_dict[rcm_name]['ssp370']['fx']}*.nc")
            lat=mask_data.lat.values
            lon=mask_data.lon.values
            abs_diff = np.abs(lat - lat_center) + np.abs(lon - lon_center)
            iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

            # Slice area around the city
            grid_number = 90 if rcm_name != 'BCCR-UCAN_eur12' else 30
            grid_number = grid_number * grid_multiplier
            half_box = grid_number // 2
            y_slice = slice(max(0, iy - half_box), iy + half_box)
            x_slice = slice(max(0, ix - half_box), ix + half_box)

            mask = mask_data['sftlf'].isel(y=y_slice, x=x_slice)
            mask_dict[rcm_name] = mask

        
    for season in seasons:
        for rcm_name in objective_rcms:
            matches_hist = data_hist[season][rcm_name]['day'].dt.floor('D') == data_hist[season][rcm_name]['1hr'].dt.floor('D')
            percentage_matches_hist[season][rcm_name]['Percentage'] = (matches_hist.sum(dim='season_year') / matches_hist.sizes['season_year']) * 100

            matches_ssp370 = data_ssp370[season][rcm_name]['day'].dt.floor('D') == data_ssp370[season][rcm_name]['1hr'].dt.floor('D')
            percentage_matches_ssp370[season][rcm_name]['Percentage'] = (matches_ssp370.sum(dim='season_year') / matches_ssp370.sizes['season_year']) * 100

    
    value_limits = [(0, 100)]
    colobar_limits = [ (0, 20, 20)]
    if custom == 1:
        custom_colorbar = [utils.custom_color_1]
        color_extra_name = '_customcolor1'
        custom_ticks = None
    elif custom == 2:
        raise NotImplementedError("No custom 2 for this plot")
        custom_colorbar = [utils.custom_color_2]
        color_extra_name = '_customcolor2'
        custom_ticks =  [utils.custom_ticks_2_day]
    elif custom == 3:
        raise NotImplementedError("No custom 3 for this plot")
        custom_colorbar = [utils.custom_color_3]
        color_extra_name = '_customcolor3'
        custom_ticks = [utils.custom_ticks_3]

    for season in seasons:
        utils.multi_map(data=percentage_matches_hist[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits,
                custom_color=custom_colorbar, cbar_limits=colobar_limits, custom_cbar_ticks=custom_ticks,
                title=f'Historical Percentage-{season} CNRM-MF vs BCCR-UCAN vs BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Percentage_Historical_date_max_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        utils.multi_map(data=percentage_matches_ssp370[season], x_map=rcm_dict, y_map=frecuency_dict, vlimits=value_limits,
                custom_color=custom_colorbar, cbar_limits=colobar_limits, custom_cbar_ticks=custom_ticks,
                title=f'SSP370 Percentage-{season}  BCCR-UCAN_BCCR-UCAN_eur12',
                fig_path=FIGURES_PATH, fig_name=f'Percentage_Future_date_max_{season}_CMIP6_{target_coords}{color_extra_name}.png',
                mask_dict = mask_dict)

        
    print("Plots terminados")
