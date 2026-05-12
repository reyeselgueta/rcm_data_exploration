import time
import os
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

rcm_dict = {rcm_name:rcm_name for rcm_name in rcm_list}
gcm_dict = {gcm_name:gcm_name for gcm_name in gcm_list}

gcm_gwl3_years = utils.gcm_gwl3_years
relative_diff = {season: {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list} for season in seasons}
relative_rx1day = {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list}
relative_prhmax = {rcm_name:{gcm_name:None for gcm_name in gcm_list} for rcm_name in rcm_list}


#prh_ref = xr.open_dataset(f'{DATA_PATH_METRICS}/hist_climatology_{metric}_prhmax_CNRM-CM5_CCLM4-8-17_1986-2005.nc')
rx1day_ref = xr.open_dataset(f'{DATA_PATH_METRICS}/hist_climatology_{metric}_rx1day_CNRM-CM5_CCLM4-8-17_1986-2005.nc')
for rcm_name in rcm_list:
    for gcm_name in gcm_list:
        # Load metrics
        path_hist_rx1day = f'{DATA_PATH_METRICS}/hist_climatology_{metric}_rx1day_{gcm_name}_{rcm_name}_1986-2005.nc'
        path_fut_rx1day = f'{DATA_PATH_METRICS}/gwl3_climatology_{metric}_rx1day_{gcm_name}_{rcm_name}_{gcm_gwl3_years[gcm_name][0]}-{gcm_gwl3_years[gcm_name][1]}.nc'
        if not Path(path_hist_rx1day).is_file() or not Path(path_fut_rx1day).is_file():
            print(f"⚠️ Missing files for {gcm_name} + {rcm_name} - RX1DAY")
            continue
        metric_hist_rx1day = xr.open_dataset(path_hist_rx1day)
        metric_fut_rx1day = xr.open_dataset(path_fut_rx1day)
        metric_hist_rx1day = metric_hist_rx1day.rename({'rx1day': metric})
        metric_fut_rx1day = metric_fut_rx1day.rename({'rx1day': metric})

        relative_rx1day_temp = (metric_fut_rx1day - metric_hist_rx1day) / metric_hist_rx1day * 100
        relative_rx1day[rcm_name][gcm_name] = relative_rx1day_temp

        path_hist_prhmax = f'{DATA_PATH_METRICS}/hist_climatology_{metric}_prhmax_{gcm_name}_{rcm_name}_1986-2005.nc'
        path_fut_prhmax = f'{DATA_PATH_METRICS}/gwl3_climatology_{metric}_prhmax_{gcm_name}_{rcm_name}_{gcm_gwl3_years[gcm_name][0]}-{gcm_gwl3_years[gcm_name][1]}.nc'
        if not os.path.exists(path_hist_prhmax) or not os.path.exists(path_fut_prhmax):
            print(f"⚠️ Missing files for {gcm_name} + {rcm_name} in prhmax")
            continue

        metric_hist_prhmax = xr.open_dataset(path_hist_prhmax)
        metric_hist_prhmax = utils.fix_latlon(metric_hist_prhmax)
        metric_hist_prhmax = metric_hist_prhmax.rename({'prhmax': metric})
        metric_hist_prhmax = utils.smart_regrid(metric_hist_prhmax, rx1day_ref, var=metric)
        
        metric_fut_prhmax = xr.open_dataset(path_fut_prhmax)
        metric_fut_prhmax = utils.fix_latlon(metric_fut_prhmax)
        metric_fut_prhmax = metric_fut_prhmax.rename({'prhmax': metric})
        metric_fut_prhmax = utils.smart_regrid(metric_fut_prhmax, rx1day_ref, var=metric)
        
        relative_prhmax_temp = (metric_fut_prhmax - metric_hist_prhmax) / metric_hist_prhmax * 100
        relative_prhmax[rcm_name][gcm_name] = relative_prhmax_temp

        for season in seasons:
            relative_diff[season][rcm_name][gcm_name] = relative_prhmax_temp.sel(season=season) - relative_rx1day_temp.sel(season=season)
        # print("Relative difference (PRHMAX - RX1DAY)")
        # print(relative_diff[rcm_name][gcm_name])

# for season in seasons:
#     utils.multi_map(data=relative_diff[season], x_map=rcm_dict, y_map=gcm_dict, vlimits=(-50, 50), var=metric,
#             color='BrBG', cbar_limits=(0, 10, 10), title=f'Relative difference {season} prhmax - rx1day ({metric}) - % (Average rx1day over 20 years, 1986-2005 as reference, and gwl3 as target.))',
#             fig_path=FIG_PATH, fig_name=f'Difference_Relatives_{metric}_{season}.png')
    

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

def plot_seasonal_matrix(data_dict, gcm_list, rcm_list,
                        var_name='yr-mean', title='None', fig_name='Nombre.png',
                        color='BrBG'):
    """
    Plotea una matriz de mapas para una season específica.
    Columnas: GCMs | Filas: RCMs
    """
    n_rows = len(rcm_list)
    n_cols = len(gcm_list)
    
    # Crear la figura con proyección Cartopy
    # Ajustamos el tamaño según la cantidad de modelos
    fig, axes = plt.subplots(
        n_rows, n_cols, 
        figsize=(n_cols * 5.2, n_rows * 3.2),
        subplot_kw={'projection': ccrs.PlateCarree()},
        sharex=True, sharey=True,
        constrained_layout=True
    )
    
    # Definir niveles para un colorbar discreto (BrBG suele ser para diferencias)
    # Puedes ajustar estos niveles según el rango de tus datos
    levels = np.linspace(-100, 100, 11) # Ejemplo: de -100% a 100% con 10 saltos
    cmap = plt.get_cmap(color, len(levels) - 1)

    for r_idx, rcm in enumerate(rcm_list):
        for g_idx, gcm in enumerate(gcm_list):
            ax = axes[r_idx, g_idx]
            
            # Obtener el dataset del diccionario anidado
            ds = data_dict[rcm][gcm]
            
            if ds is not None:
                # Plot de los datos
                im = ds[var_name].plot(
                    ax=ax, 
                    transform=ccrs.PlateCarree(),
                    levels=levels,
                    cmap=cmap,
                    add_colorbar=False, # Quitamos colorbars individuales
                    add_labels=False    # Quitamos etiquetas de ejes internas
                )
                
                # Añadir detalles geográficos
                ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
                ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
                ax.set_extent([ds.lon.min(), ds.lon.max(), ds.lat.min(), ds.lat.max()])
            
            # Títulos solo en la primera fila (GCMs)
            if r_idx == 0:
                ax.set_title(f"GCM: {gcm}", fontweight='bold', fontsize=16)
            else:
                ax.set_title("")
                
            # Etiquetas solo en la primera columna (RCMs)
            if g_idx == 0:
                ax.text(-0.2, 0.5, f"RCM: {rcm}", transform=ax.transAxes, 
                        rotation=90, va='center', ha='right', fontweight='bold', fontsize=16)

    # Añadir el Colorbar común al final de todos los subplots
    cbar = fig.colorbar(
        im, ax=axes, orientation='vertical', 
        shrink=1.0, pad=0.03, aspect=30, extend='both'
    )
    cbar.set_label(f"Relative Difference (%)", fontsize=20, fontweight='bold')
    cbar.ax.tick_params(labelsize=16) # Tamaño números colorbar

    plt.suptitle(f"{title}", fontsize=24, fontweight='bold', y=1.05)
    plt.savefig(f'{FIG_PATH}/{fig_name}', bbox_inches='tight')
    plt.close()

# --- Ejemplo de uso ---
for season in seasons:
    plot_seasonal_matrix(relative_diff[season], gcm_list, rcm_list,
                        color = 'BrBG',
                        title = f"Relative diff - metric:{metric} - season_{season}",
                        fig_name =f"Difference_Relatives_{metric}_{season}.png" 
                        )
