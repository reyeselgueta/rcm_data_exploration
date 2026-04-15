import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
from typing import List, Tuple, Optional

class TimeSeriesPlot:
    """
    This class provides a structured way to create time series plots for multiple metrics, RCMs, and frequencies. It allows for consistent styling and easy addition of subplots, as well as a global legend that encompasses all the plotted data.

    """
    def __init__(self, n_rows, n_cols, figsize_scale=4):
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.fig, self.axes = plt.subplots(
            self.n_rows,
            self.n_cols,
            figsize=(self.n_cols * figsize_scale, self.n_rows * figsize_scale),
        )
        # Diccionario global de handles y labels de la leyenda
        self.legend_handles = {}
    
    def plot_subplot(
        self,
        row_idx,
        col_idx,
        data_metric,  # [{'data_list': [...], 'label': label, 'color': color}]
        time_var,
        ylims=None,
        title="",
        y_label="",
        x_label="",
    ):
        if self.n_cols == 1:
            ax = self.axes[row_idx]
        else:
            ax = self.axes[row_idx, col_idx]

        for data_dict in data_metric:
            label = data_dict['label']
            color = data_dict['color']
            for data_element in data_dict['data_list']:
                # if isinstance(data_element, xr.Dataset):
                #     var_name = list(data_element.data_vars)[0]
                #     data_element = data_element[var_name]
                if isinstance(data_element, xr.Dataset):
                    data_element = data_element.to_dataarray().squeeze()
                #ax.plot(data_element[time_var], data_element, label=label, color=color, lw=1.5)
                # data_element['time'] = xr.DataArray(
                #     np.array(data_element['time'].values, dtype='datetime64[ns]'),
                #     dims='time'
                # )
                ax.plot(
                    data_element[time_var],
                    data_element,
                    label=label,
                    color=color,
                    lw=1.5,
                    marker=None
                )
            # Guardar un handle para la leyenda global
            # if label not in self.legend_handles:
            #     line = plt.Line2D([0], [0], color=color, lw=2)
            #     self.legend_handles[label] = line

        ax.set_title(title, fontsize=9)
        if ylims is not None:
            ax.set_ylim(ylims)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.set_ylabel(y_label, fontsize=9)
        ax.set_xlabel(x_label, fontsize=9)

    def finalize(self, savepath, title=""):
        # Crear handles de la leyenda para todos los colores definidos en COLORS
        handles = [plt.Line2D([0], [0], color=color, lw=2) for color in colors.values()]
        labels = list(colors.keys())

        self.fig.legend(
            handles=handles,
            labels=labels,
            loc="upper right",
            ncol=int(np.ceil(len(labels) / 5)),
            fontsize=8
        )
        self.fig.suptitle(title, fontsize=12)
        self.fig.tight_layout(rect=[0, 0, 1, 0.97])
        self.fig.savefig(savepath)
        #plt.close(self.fig)





colors = {
    'CNRM-MF': 'blue',
    'BCCR-UCAN': 'orange',
    'BCCR-UCAN_eur12': 'green'
    }