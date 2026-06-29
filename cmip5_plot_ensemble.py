import time
import xarray as xr
import numpy as np
from pathlib import Path
import sys
import utils_explore as utils


DATA_PATH_METRICS = '/nfs/home/gmeteo/reyess/rcm_exploration/data/metrics/'
FIG_PATH = '/nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration/example_figures/ensemble/'
metric = sys.argv[1] # mon-mean, yr-mean, max-mean
area = sys.argv[2] #Valencia or Iberia
# metric = "yr-mean" # mon-mean, yr-mean, max-mean
# area = "Iberia"
print(f"Doing metric: {metric} - area: {area}")

lat_target, lon_target = utils.get_slice_coords(area)


rcm_list = ['CCLM4-8-17', 'COSMO-crCLIM-v1-1', 'HIRHAM5',
            'HadREM3-GA7-05', 'RACMO22E']
gcm_list = ['CNRM-CM5', 'EC-EARTH', 'HadGEM2-ES', 'MPI-ESM-LR', 'NorESM1-M']
seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']
row_metrics = ['Mean', 'P80', 'P20']


rcm_dict = {rcm_name:rcm_name for rcm_name in rcm_list}
gcm_dict = {gcm_name:gcm_name for gcm_name in gcm_list}

gcm_gwl3_years = utils.gcm_gwl3_years
prhmax_data = {row_metric: {season: None for season in seasons} for row_metric in row_metrics}
rx1day_data = {row_metric: {season: None for season in seasons} for row_metric in row_metrics}

# Load metrics
for row_metric in row_metrics:
    ds_rx1day = xr.open_dataset(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_{row_metric.lower()}_{area}_{metric}_rx1day.nc')
    ds_prhmax = xr.open_dataset(f'{DATA_PATH_METRICS}/cmip5ensemble_relative-ensemble_{row_metric.lower()}_{area}_{metric}_prhmax.nc')

    for season in seasons:
        rx1day_data[row_metric][season] = ds_rx1day.sel(season=season)
        prhmax_data[row_metric][season] = ds_prhmax.sel(season=season)


# utils.multi_map(data=rx1day_data, x_map=row_metrics, y_map=seasons, vlimits=[(0, 30), (0, 30), (0, 30)], var='rx1day',
#         color=['BrBG', 'BrBG', 'BrBG'], cbar_limits=[(0, 10, 10), (0, 10, 10), (0, 10, 10)], title=f'Metrics rx1day {metric}',
#         fig_path=FIG_PATH, fig_name=f'Metric_Rx1day_{metric}.png')

# utils.multi_map(data=prhmax_data, x_map=row_metrics, y_map=seasons, vlimits=[(0, 30), (0, 30), (0, 30)], var='prhmax',
#         color=['BrBG', 'BrBG', 'BrBG'], cbar_limits=[(0, 10, 10), (0, 10, 10), (0, 10, 10)], title=f'Metrics prhmax {metric}',
#         fig_path=FIG_PATH, fig_name=f'Metric_prhmax_{metric}.png')

import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

def plot_full_matrix(data_dict, seasons_list, metrics_order=['Mean', 'P20', 'P80'],
                     save_path='matriz_precipitacion.png', 
                     title='None', y_labels=None, x_labels=None,
                     color='Blues', vlimits=None):
    """
    Plotea una matriz completa de mapas. Corrige el granulado detectando
    matrices de coordenadas 2D (latitude/longitude) o vectores 1D (lat/lon).
    """
    n_rows = len(metrics_order)
    n_cols = len(seasons_list)
    proyeccion = ccrs.PlateCarree() 
    
    primera_metrica = metrics_order[0]
    primera_season = seasons_list[0]
    sample_ds = data_dict.get(primera_metrica, {}).get(primera_season, None)
    
    if sample_ds is not None:
        data_vars = list(sample_ds.data_vars)
        var_name = data_vars[0] if data_vars else 'rx1day'
    else:
        var_name = 'rx1day'

    colorbars_por_fila = isinstance(color, list)

    if y_labels is None:
        y_labels = metrics_order
    if x_labels is None:
        x_labels = seasons_list

    fig, axes = plt.subplots(
        n_rows, n_cols, 
        figsize=(n_cols * 3.5, n_rows * 3.4),
        subplot_kw={'projection': proyeccion},
        sharex=True, sharey=True,
        constrained_layout=True
    )

    im_por_fila = {}
    im_general = None

    for r_idx, metric in enumerate(metrics_order):
        row_color = color[r_idx] if colorbars_por_fila else color
        
        if vlimits is None:
            seasons_to_check = seasons_list
            metrics_to_check = [metric] if colorbars_por_fila else metrics_order
            all_values = [data_dict[m][s][var_name].values 
                          for m in metrics_to_check for s in seasons_to_check 
                          if m in data_dict and s in data_dict[m] and var_name in data_dict[m][s]]
            row_vmin = np.nanmin(all_values) if all_values else 0
            row_vmax = np.nanmax(all_values) if all_values else 100
        else:
            row_vmin, row_vmax = vlimits[r_idx] if isinstance(vlimits, list) else vlimits

        levels = np.linspace(row_vmin, row_vmax, 11)
        cmap = plt.get_cmap(row_color, len(levels) - 1)

        for c_idx, season in enumerate(seasons_list):
            ax = axes[r_idx, c_idx]
            ds = data_dict.get(metric, {}).get(season, None)
            
            if ds is not None:
                # --- DETECTOR INTELIGENTE DE COORDENADAS 2D ---
                # Si el dataset tiene las matrices completas de 'longitude' y 'latitude'
                if 'longitude' in ds.coords and 'latitude' in ds.coords:
                    x_dim, y_dim = 'longitude', 'latitude'
                    plot_transform = ccrs.PlateCarree() # Al usar lat/lon reales, el transform es PlateCarree
                elif 'lon' in ds.coords and 'lat' in ds.coords:
                    x_dim, y_dim = 'lon', 'lat'
                    plot_transform = ccrs.PlateCarree()
                else:
                    # Caso de respaldo por si vinieran rlon/rlat como vectores 1D puros
                    x_dim, y_dim = 'rlon', 'rlat'
                    if 'rotated_pole' in ds.variables:
                        plot_transform = ccrs.RotatedPole(
                            pole_longitude=ds['rotated_pole'].attrs.get('grid_north_pole_longitude', 0),
                            pole_latitude=ds['rotated_pole'].attrs.get('grid_north_pole_latitude', 90)
                        )
                    else:
                        plot_transform = ccrs.RotatedPole(pole_longitude=-162.0, pole_latitude=39.25)

                # Pintar los datos
                im = ds[var_name].plot(
                    ax=ax, transform=plot_transform,
                    x=x_dim, y=y_dim,
                    levels=levels, cmap=cmap,
                    add_colorbar=False, add_labels=False
                )
                im_por_fila[r_idx] = im
                im_general = im
                
                ax.add_feature(cfeature.COASTLINE, linewidth=0.6, edgecolor='black')
                ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
                
                # Para evitar desajustes en el zoom debido a los NaN de los bordes,
                # calculamos el extent usando los valores reales y seguros de la proyección de dibujo
                ax.set_extent([ds[x_dim].min(), ds[x_dim].max(), 
                               ds[y_dim].min(), ds[y_dim].max()], crs=plot_transform)
            
            if r_idx == 0:
                ax.set_title(x_labels[c_idx], fontweight='bold', fontsize=14)
            else:
                ax.set_title("")
                
            if c_idx == 0:
                ax.text(-0.25, 0.5, y_labels[r_idx], transform=ax.transAxes, 
                        rotation=90, va='center', ha='right', fontweight='bold', fontsize=14)

    if colorbars_por_fila:
        for r_idx in range(n_rows):
            if r_idx in im_por_fila:
                cbar = fig.colorbar(
                    im_por_fila[r_idx], ax=axes[r_idx, :], 
                    orientation='horizontal', shrink=0.6, pad=0.04, extend='both'
                )
                cbar.ax.tick_params(labelsize=10)
    else:
        if im_general is not None:
            cbar = fig.colorbar(im_general, ax=axes, orientation='vertical', shrink=0.7, pad=0.02, extend='both')
            cbar.set_label(f"Value ({var_name})", fontsize=14, fontweight='bold')

    plt.suptitle(title, fontsize=18, fontweight='bold', y=1.03)
    
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, bbox_inches='tight', dpi=150)
    plt.close()
# ==========================================
# --- Ejemplo de ejecución único ---
# ==========================================

seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']

# Nombres más amigables en español para los ejes
eje_x_estaciones = ['Anual', 'Invierno (DJF)', 'Primavera (MAM)', 'Verano (JJA)', 'Otoño (SON)']
eje_y_metricas = ['Media', 'Percentil 20 (P20)', 'Percentile 80 (P80)']

mean_limits = (0, 40) if metric != 'max-mean' else (0,90)
p20_limits = (0, 8)
p80_limits = (0, 80) if metric != 'max-mean' else (0,180)
vlimits = [mean_limits, p20_limits, p80_limits]  # Límite de valores para cada métrica
# Llamada para rx1day (Detectará la variable sola y aplicará sus vlimits correspondientes)
plot_full_matrix(
    data_dict=rx1day_data, 
    seasons_list=seasons, 
    metrics_order=['Mean', 'P20', 'P80'],
    save_path=f"{FIG_PATH}/rx1day_metrics_{metric}_{area}.png",
    title=f"Rx1day metrics ({metric})",
    x_labels=eje_x_estaciones,
    y_labels=eje_y_metricas,
    color=['YlGnBu', 'YlGnBu', 'YlGnBu'],  
    vlimits=vlimits  
)

# Llamada para prhmax (Detectará prhmax sola, adaptará los ejes rlat/rlon y usará estos vlimits)
plot_full_matrix(
    data_dict=prhmax_data, 
    seasons_list=seasons, 
    metrics_order=['Mean', 'P20', 'P80'],
    save_path=f"{FIG_PATH}/prhmax_metrics_{metric}_{area}.png",
    title=f"Prhmax metrics ({metric})",
    x_labels=eje_x_estaciones,
    y_labels=eje_y_metricas,
    color=['YlGnBu', 'YlGnBu', 'YlGnBu'],  
    vlimits=vlimits 
)
