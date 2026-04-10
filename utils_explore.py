import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
from typing import List, Tuple, Optional

def multi_map(
    data: xr.Dataset,
    x_map: dict,
    y_map: dict,
    fig_path: str,
    fig_name: str,
    vlimits: Optional[List[Tuple[float, float]]] = None,
    color: List[str] = ['RdBu_r'],
    custom_color=None,
    cbar_limits: List[Tuple[int, int, int]] = [(0, 10, 10)],
    custom_cbar_ticks=None,
    var: Optional[str] = None,
    lon: str = 'lon',
    lat: str = 'lat',
    title: str = '',
    mask_dict = None
):
    """
    Create a multi-panel map plot.

    Parameters:
    - data: The dataset to plot.
    - x_map: A dictionary mapping x-axis labels to their corresponding values.
    - y_map: A dictionary mapping y-axis labels to their corresponding values.
    - fig_path: The path where the figure will be saved.
    - fig_name: The name of the figure file.
    - vlimits: A list of tuples specifying the value limits for each subplot.
    - color: A list of colormap names.
    - custom_color: A list of custom colormap colors.
    - cbar_limits: A list of tuples specifying the colorbar limits for each subplot.
    - custom_cbar_ticks: A list of custom colorbar ticks.
    - var: The variable to plot.
    - lon: The longitude coordinate name.
    - lat: The latitude coordinate name.
    - title: The title for the figure.
    - mask_dict: A dictionary containing mask data.

    Returns:
    None
    """

    n_rows = len(y_map.keys())
    n_cols = len(x_map.keys())

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4*n_cols, 3*n_rows),
        sharex=False, sharey=False
    )

    # asegurar que axes es 2D
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes[np.newaxis, :]
    elif n_cols == 1:
        axes = axes[:, np.newaxis]

    row_mappables = []
    
    for j, temporal_res in enumerate(y_map.keys()):
        is_over_max = False
        is_under_min = False
        row_maps = []

        for i, rcm_name in enumerate(x_map.keys()):
            ax = axes[j, i]

            if custom_cbar_ticks is not None:
                cmap = mcolors.ListedColormap(custom_color[j])
                norm = mcolors.BoundaryNorm(boundaries=custom_cbar_ticks[j], ncolors=len(custom_color[j]))
            else:
                continuousCMAP = plt.get_cmap(color[j]) if custom_color is None else ListedColormap(custom_color[j])
                cmap = ListedColormap(
                    continuousCMAP(np.linspace(*(0, 1, cbar_limits[j][2]))[cbar_limits[j][0]:cbar_limits[j][1]])
                )
                norm = None

            if j == 0 and x_map is not None:
                ax.set_title(f'{x_map[rcm_name]}', fontsize=12, pad=10)
            if i == 0 and y_map is not None:
                ax.text(-0.20, 0.55, f'{y_map[temporal_res]}', va='bottom', ha='center',
                        rotation='vertical', rotation_mode='anchor',
                        transform=ax.transAxes, fontsize=12)

            dataToPlot = data[rcm_name][temporal_res][var] if var is not None else data[rcm_name][temporal_res]

            # ploteo principal, sin transform
            im = ax.pcolormesh(
                dataToPlot.coords[lon].values,
                dataToPlot.coords[lat].values,
                dataToPlot,
                cmap=cmap,
                norm=norm,
                vmin=None if norm else vlimits[j][0],
                vmax=None if norm else vlimits[j][1],
                shading='auto'
            )
            row_maps.append(im)

            # dibujar contorno de máscara
            if mask_dict is not None:
                ax.contour(
                    mask_dict[rcm_name][lon], mask_dict[rcm_name][lat], mask_dict[rcm_name].values,
                    levels=[0.5],
                    colors='black',
                    linewidths=1.0
                )

            # valor medio en esquina
            median_val = float(dataToPlot.median().values)
            ax.text(
                0.98, 0.10, f"median: {median_val:.2f}",
                transform=ax.transAxes,
                fontsize=6, color='black',
                ha='right', va='bottom',
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
            )
            p95 = float(np.nanpercentile(dataToPlot.values, 95))

            ax.text(
                0.98, 0.18, f"p95: {p95:.2f}",
                transform=ax.transAxes,
                fontsize=6, color='black',
                ha='right', va='bottom',
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
            )
            p5 = float(np.nanpercentile(dataToPlot.values, 5))
            ax.text(
                0.98, 0.02, f"p5: {p5:.2f}",
                transform=ax.transAxes,
                fontsize=6, color='black',
                ha='right', va='bottom',
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
            )
            if dataToPlot.max().values > vlimits[j][1]:
                is_over_max = True
            if dataToPlot.min().values < vlimits[j][0]:
                is_under_min = True

        row_mappables.append(row_maps[-1])

        # colorbar por fila
        if is_over_max and is_under_min:
            extend = 'both'
        elif is_under_min:
            extend = 'min'
        elif is_over_max:
            extend = 'max'
        else:
            extend = 'neither'
        cbar = fig.colorbar(
            row_maps[-1],
            ax=axes[j, :],
            orientation='vertical',
            fraction=0.06,
            pad=0.06,
            extend=extend
        )
        if custom_cbar_ticks is not None:
            cbar.set_ticks(custom_cbar_ticks[j])
            cbar.set_ticklabels([str(v) for v in custom_cbar_ticks[j]])
        else:
            ticks = np.linspace(vlimits[j][0], vlimits[j][1], int(np.floor(cbar_limits[j][1] - cbar_limits[j][0]))+1)
            cbar.set_ticks(ticks)
        cbar.ax.tick_params(labelsize=10)

    plt.suptitle(title, fontsize=12, y=1.01)
    plt.savefig(f'{fig_path}/{fig_name}', bbox_inches='tight')
    plt.close()


def plot_data(
    datain,
    mask_dict,
    fig_path,
    fig_name,
    title="",
    varname="pr",
    lon_city=-0.375, lat_city=39.467,
    n_grid=30,
    vmin=[0, 0],
    vmax=[1200, 50],
    cmap=['viridis', 'BrBG'],
    name="max"
):
    """Plot a single map of data with mask and city marker."""

    temporalities = list(next(iter(datain.values())).keys())  # ['day','1hr']
    rcms = list(datain.keys())  # ['rcm1','rcm2']

    n_rows = len(temporalities)
    n_cols = len(rcms)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(5*n_cols, 4*n_rows),
        sharex=True, sharey=True,
        squeeze=False
    )

    for j, temp in enumerate(temporalities):   # filas
        row_images = []  # para guardar handles de pcolormesh o plot
        for i, rcm in enumerate(rcms):         # columnas
            ax = axes[j, i]

            data = datain[rcm][temp]#[varname]

            im = data.plot(
                ax=ax,
                x="lon", y="lat",
                cmap=cmap[j],
                vmin=vmin[j], vmax=vmax[j],
                add_colorbar=False  # colorbar manual después
            )
            row_images.append(im)

            # Contorno máscara
            mask_dict[rcm].plot.contour(
                ax=ax,
                x="lon", y="lat",
                levels=[50],
                colors="black",
                linewidths=1.0,
                zorder=10
            )

            # Etiquetas de columnas y filas
            if j == 0:
                ax.set_title(rcm, fontsize=12)
            if i == 0:
                ax.set_ylabel(temp, fontsize=12)

        # Colorbar para la fila entera (a la derecha)
        cbar = fig.colorbar(
            row_images[0],
            ax=axes[j, :],  # toda la fila
            orientation="vertical",
            fraction=0.05,
            pad=0.08
        )
        cbar.set_label(varname, fontsize=10)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout(rect=[0, 0, 0.95, 0.95])

    os.makedirs(fig_path, exist_ok=True)
    outfile = f"{fig_path}/{fig_name}"
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Figure saved to: {outfile}")



def multi_map_ccrs(data:xr.Dataset, x_map:dict, y_map:dict, fig_path:str, fig_name:str,
            vlimits:tuple=None, color:str='RdBu_r', custom_color=None, cbar_limits:tuple=(0, 10, 10),
            var=None, lon:str='lon', lat:str='lat', title:str=''):

    n_rows = len(y_map.keys())
    n_cols = len(x_map.keys())

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), sharex=False, sharey=False, subplot_kw={'projection': ccrs.PlateCarree()})
    continuousCMAP = plt.get_cmap(color) if custom_color==None else ListedColormap(custom_color)
    discreteCMAP = ListedColormap(continuousCMAP(np.linspace(*(0, 1, cbar_limits[2]))[cbar_limits[0]:cbar_limits[1]]))
    is_over_max = False
    is_under_min = False
    for i, rcm_name in enumerate(x_map.keys()):
        for j, gcm_name in enumerate(y_map.keys()):
            print(f"{rcm_name} - GCM: {gcm_name}")
            if n_cols == 1 and n_rows == 1:
                ax = axes
            elif n_cols == 1:
                ax = axes[j]
            elif n_rows == 1:
                ax = axes[i]
            else:
                ax = axes[j, i]


            if j == 0 and x_map is not None:
                    ax.set_title(f'{x_map[rcm_name]}', fontsize=12, pad=10)
            if i == 0 and y_map is not None:
                ax.text(-0.07, 0.55, f'{y_map[gcm_name]}', va='bottom', ha='center',
                    rotation='vertical', rotation_mode='anchor',
                    transform=ax.transAxes, fontsize=12)
            ax.coastlines(resolution='10m')    

            if data[rcm_name][gcm_name] is None:
                ax.set_axis_off()
                continue
            dataToPlot = data[rcm_name][gcm_name][var] if var!=None else data[rcm_name][gcm_name]#.isel(member=0, drop=True)#[var]
            if 'newlon' in dataToPlot.coords:
                current_lat = 'newlat'
                current_lon = 'newlon'
            else:
                current_lat = lat
                current_lon = lon
            dataAggregated = dataToPlot.mean(dim=list(dataToPlot.dims))
            im = ax.pcolormesh(dataToPlot.coords[current_lon].values, dataToPlot.coords[current_lat].values,
                                dataToPlot,
                                #transform=ccrs.PlateCarree(),
                                cmap=discreteCMAP,
                                vmin=vlimits[0], vmax=vlimits[1])
            
            # Agregar el valor medio (dataAggregated) en la esquina inferior derecha
            mean_val = float(dataAggregated.values)
            ax.text(
                0.98, 0.02, f"{mean_val:.2f}",  
                transform=ax.transAxes,       
                fontsize=12, color='black',
                ha='right', va='bottom',     
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')  
            )
            if dataToPlot.max().compute() > vlimits[1]:
                is_over_max = True
            if dataToPlot.min().compute() < vlimits[0]:
                is_under_min = True

    cax = fig.add_axes([0.925, 0.125, 0.02, 0.775])
    if is_over_max and is_under_min:
        extend = 'both'
    elif is_under_min:
        extend = 'min'
    elif is_over_max:
        extend = 'max'
    else:
        extend = 'neither'
    cbar = plt.colorbar(im, cax, orientation='vertical', spacing='uniform', extend=extend)
    ticks = np.linspace(vlimits[0], vlimits[1], int(np.floor(cbar_limits[1] - cbar_limits[0]))+1)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=12)
    plt.suptitle(title, fontsize=12)
    plt.subplots_adjust(top=0.95, bottom=0.05, wspace=0.002, hspace=0.002)
    plt.savefig(f'{fig_path}/{fig_name}', bbox_inches='tight')
    plt.close()


def get_season_data(data, season, year_to_drop=2015, chunk=False):
    
    if chunk == True:
        data = data.chunk({"time": 1000})

    if season != 'Annual':
        if season == 'DJF':
            # data.coords["season_year"] = data["time"].dt.year
            # data["season_year"] = xr.where(data["time"].dt.month == 12,
            #                             data["season_year"] + 1,
            #                             data["season_year"])
            season_year = data["time"].dt.year
            season_year = xr.where(data["time"].dt.month == 12,
                                season_year + 1,
                                season_year)
            data = data.assign_coords(season_year=season_year)
            data_yearly = data.where((data["time.season"] == season) & (data.season_year != year_to_drop)
                                                        , drop=True).groupby("season_year")
        else:
            data_yearly = data.where(data["time.season"] == season, drop=True).groupby("time.year")
    else:
        data_yearly = data.groupby("time.year")
        
    return data_yearly

def get_statistic(data, statistic, time_multiplier=1, season = None, spatial = True):
    """
    Calculate the specified statistic from the data.
    Parameters:
    - data: xarray Dataset or DataArray containing the data to analyze. For spatial is yearly data grouped by season_year or time.year, for aggregated is the original data.
    - statistic: string specifying the statistic to calculate. Can be 'mean', 'max', 'P95', 'P995', or 'date_max'.
    - time_multiplier: numeric value to multiply the result by, used for unit conversion (e.g., from kg/m²/s to mm/day).
    - season: string specifying the season (e.g., 'DJF', 'MAM', 'JJA', 'SON', 'Annual'), used to determine the time dimension for spatial data.
    - spatial: boolean indicating whether the data is spatial (True) or aggregated (False). For spatial data, the function calculates the statistic across time and then averages or finds the max
    """
    if spatial:
        time_dim = 'year' if season != 'DJF' else 'season_year'

        if statistic == 'mean-max':
            result = data.max(dim="time").mean(dim=time_dim)

        elif statistic == 'max-max':
            result = data.max(dim="time").max(dim=time_dim)

        elif statistic == 'mean-P95':
            result = data.quantile(0.95, dim="time").mean(dim=time_dim)

        elif statistic == 'mean-P995':
            result = data.quantile(0.995, dim="time").mean(dim=time_dim)

        elif statistic == 'date_max':
            idx = data.idxmax(dim="time")
            result = xr.decode_cf(idx.to_dataset(name="date_max"))

        elif statistic == 'mean':
            result = data.mean(dim="time").mean(dim=time_dim)

        else:
            raise ValueError(f"Unknown statistic: {statistic}")

    elif spatial == False:

        if statistic == 'mean':
            result = data.mean(dim=("y", "x"))

        elif statistic == 'max':
            result = data.max(dim=("y", "x"))

        elif statistic == 'date_max':
            # 🔥 primero colapsar espacio
            spatial_max = data.max(dim=("y", "x"))

            # luego encontrar el tiempo del máximo
            idx = spatial_max.idxmax(dim="time")

            result = xr.decode_cf(idx.to_dataset(name="date_max"))

        else:
            raise ValueError(f"Unknown statistic: {statistic}")

    return result * time_multiplier

def get_coords(target_coords):
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
    
    return lat_center, lon_center, folder_name, grid_multiplier

gcm_gwl3_years = {'CNRM-CM5': ('2058', '2077'), 'CanESM2': ('2040', '2059'),
            'EC-EARTH': ('2051', '2070'), 'HadGEM2-ES': ('2045', '2064'), 'IPSL-CM5A-MR': ('2041', '2060'),
            'MIROC5': ('2063', '2082'), 'MPI-ESM-LR': ('2052', '2071'), 'NorESM1-M': ('2063', '2082'),
            'NorESM2-MM': ('2081', '2100'), 'CNRM-ESM2-1': ('2072', '2091'), 'EC-Earth3-Veg': ('2052', '2071')} #Confirmar estos ultimos 3

custom_color_1 = [
    "#ffffff", "#e6f7ff", "#cceeff", "#b3e6ff", "#99ddff", "#80d4ff",
    "#66ccff", "#4dc3ff", "#33bbff", "#1ab2ff", "#00aaff",
    "#00bfb3", "#00cc88", "#66d96d", "#c2e255", "#ffe945",
    "#ffd032", "#ffac24", "#ff7b22", "#d73027"
]

custom_color_2 = [
    "#F7F7F7", "#CCE5F6", "#99CCFF","#66FFCC","#00FF66", "#FFFF00",
    "#FFCC00","#FF9933","#FF6600","#FF0000","#CC0000","#990099","#660066"
]

custom_color_3 = ['#D3D3D3', # Grey
    '#fff7bc','#fee391','#fec44f', # Yellow group (3 colors, decreasing lightness) - inspired by YlOrBr
    '#bae4b3','#74c476','#31a354', # Green group (3 colors, decreasing lightness) - inspired by YlGn
    '#9ECAE1', '#4292C6', '#08519C', # Blue group (3 colors, decreasing lightness) - inspired by Blues
    '#9e9ac8', '#756bb1','#54278f' # Violet group (3 colors, decreasing lightness) - inspired by Purples
    ]

custom_ticks_2_1hr = np.array([0, 0.1, 1, 3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
custom_ticks_2_day = np.array([0, 2.4, 24, 72, 120, 240, 360, 480, 600, 720, 840, 960, 1080, 1200])
custom_ticks_flux = np.array([0, 0.0002, 0.0005, 0.0010, 0.0015, 0.003, 0.0045, 0.006, 0.0075, 0.009, 0.0105, 0.012, 0.0135, 0.015])
custom_ticks_not_used = np.array([0, 0.1, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30])
custom_ticks_3 = np.array([0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500])

raw_filepath = [
'Fill with URL to data'
]

data_rcm_paths = {
    'RegCM5-0': 
    {'historical': 
        {'pr':
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/VAL-3/ICTP/EC-Earth3-Veg/historical/r1i1p1f1/RegCM5-0/v1-r1/day/pr/v20250725/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/VAL-3/ICTP/EC-Earth3-Veg/historical/r1i1p1f1/RegCM5-0/v1-r1/fx/sftlf/v20250725/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/VAL-3/ICTP/EC-Earth3-Veg/historical/r1i1p1f1/RegCM5-0/v1-r1/1hr/pr/v20250725/'
            }, 
        'tas':
            None
        },
    'ssp370': None,
    'gcm': 'EC-Earth3-Veg'},
    'BCCR-UCAN': 
    {'historical': 
        {'pr':
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/day/pr/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/1hr/pr/v20240710/'
            }, 
        'tas':
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/day/tas/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/1hr/tas/v20240710/'
            }
        },
    'ssp370': 
        {'pr':
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/day/pr/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/1hr/pr/v20240710/'
            }, 
        'tas':
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/day/tas/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/1hr/tas/v20240710/'
            }
        },
    'gcm': 'NorESM2-MM'},
    'CNRM-MF': 
    {'historical': 
        {'pr':
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/CNRM-MF/CNRM-ESM2-1/historical/r1i1p1f2/CNRM-AROME46t1/v1-r1/day/pr/v20250731/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/CNRM-MF/CNRM-ESM2-1/historical/r1i1p1f2/CNRM-AROME46t1/v1-r1/fx/sftlf/v20250731/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/CNRM-MF/CNRM-ESM2-1/historical/r1i1p1f2/CNRM-AROME46t1/v1-r1/1hr/pr/v20250731/'
            },
        'tas':
            None
        },
    'ssp370': 
        {'pr':
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/CNRM-MF/CNRM-ESM2-1/ssp370/r1i1p1f2/CNRM-AROME46t1/v1-r1/day/pr/v20250731/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/CNRM-MF/CNRM-ESM2-1/ssp370/r1i1p1f2/CNRM-AROME46t1/v1-r1/fx/sftlf/v20250731/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/ALPX-3/CNRM-MF/CNRM-ESM2-1/ssp370/r1i1p1f2/CNRM-AROME46t1/v1-r1/1hr/pr/v20250731/'
            },
        'tas':
            None
        },
    'gcm': 'CNRM-ESM2-1'},
    'BCCR-UCAN_eur12': 
    {
    'historical': 
        {'pr': 
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/day/pr/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/1hr/pr/v20240710/'
            },
        'tas':        
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/day/tas/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/historical/r1i1p1f1/WRF451I/v1-r1/1hr/tas/v20240710/'
            }
        },
    'ssp370':
        {'pr': 
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/day/pr/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/1hr/pr/v20240710/'
            },
        'tas':        
            {
            'day': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/day/tas/v20240710/',
            'fx':'/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/fx/sftlf/v20240710/',
            '1hr': '/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/EUR-12/BCCR-UCAN/NorESM2-MM/ssp370/r1i1p1f1/WRF451I/v1-r1/1hr/tas/v20240710/'
            }
        },
    'gcm': 'NorESM2-MM'
    }}


import psutil, os

def mem_status(msg=""):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024**2
    print(f"[{msg}] Memoria usada: {mem:.1f} MB")


colobar_limits = {
            'climatology': [(0, 20, 20), (0, 20, 20)],
            'delta': [(0, 20, 20), (0, 20, 20)],
            'relative': [(0, 20, 20), (0, 20, 20)]
        }
value_limits = {
            'climatology': [(0,540), (0, 180)],
            'delta': [(-150,150), (-50, 50)],
            'relative': [(-100, 100), [-100, 100]],
            'relative_per_degree': [(-40, 40), [-40, 40]]
        }
rcm_dict = {'CNRM-MF': '3 Kilometers', 'BCCR-UCAN': '3 Kilometers', 'BCCR-UCAN_eur12': '12 Kilometers'}
mask_dict = {'CNRM-MF': None, 'BCCR-UCAN': None, 'BCCR-UCAN_eur12': None}
frecuency_dict = {'day': 'Daily', '1hr': 'Hourly'}