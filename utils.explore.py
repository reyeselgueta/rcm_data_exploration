import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs



def multi_map(data:xr.Dataset, x_map:dict, y_map:dict, fig_path:str, fig_name:str,
            vlimits:tuple=None, color:str='RdBu_r', cbar_limits:tuple=(0, 10, 10),
            var=None, lon:str='lon', lat:str='lat', title:str=''):

    n_rows = len(y_map.keys())
    n_cols = len(x_map.keys())

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), sharex=False, sharey=False, subplot_kw={'projection': ccrs.PlateCarree()})
    continuousCMAP = plt.get_cmap(color)
    discreteCMAP = ListedColormap(continuousCMAP(np.linspace(*(0, 1, cbar_limits[2]))[cbar_limits[0]:cbar_limits[1]]))

    for i, rcm_name in enumerate(x_map.keys()):
        for j, gcm_name in enumerate(y_map.keys()):
            print(f"{rcm_name} - GCM: {gcm_name}")
            ax = axes[j, i]

            if j == 0 and x_map is not None:
                    ax.set_title(f'{x_map[rcm_name]}', fontsize=16)
            if i == 0 and y_map is not None:
                ax.text(-0.07, 0.55, f'{y_map[gcm_name]}', va='bottom', ha='center',
                    rotation='vertical', rotation_mode='anchor',
                    transform=ax.transAxes, fontsize=16)
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
                                transform=ccrs.PlateCarree(),
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

    cax = fig.add_axes([0.925, 0.125, 0.02, 0.775])
    cbar = plt.colorbar(im, cax, orientation='vertical', spacing='uniform')
    ticks = np.linspace(vlimits[0], vlimits[1], int(np.floor(cbar_limits[1] - cbar_limits[0]))+1)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=16)
    plt.suptitle(title, fontsize=18)
    plt.subplots_adjust(top=0.95, bottom=0.05, wspace=0.002, hspace=0.002)
    plt.savefig(f'{fig_path}/{fig_name}', bbox_inches='tight')
    plt.close()




raw_filepath = [
'Fill with URL to data'
]
