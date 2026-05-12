import time
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import os

FIGURES_PATH = '/nfs/home/gmeteo/reyess/rcm_exploration/rcm_data_exploration/example_figures/'
DATA_PATH_METRICS = '/nfs/home/gmeteo/reyess/rcm_exploration/data/metrics/'

datasets = {
    'EOBS': 'daily',
    'ROCIO': 'daily',
    'ERA5': 'hourly',
    'ERA5-Land': 'daily',
    'CERRA': 'daily',
    'CERRA-Land': 'daily'
}

grid_order = [
    ['EOBS', 'ROCIO'],
    ['ERA5', 'ERA5-Land'],
    ['CERRA', 'CERRA-Land']
]

lon = (-1.75, 1.05)
lat = (38.45, 40.55)

time_start = time.time()
# -----------------------------
# FIGURE + GRID SPEC (KEY PART)
# -----------------------------
fig = plt.figure(figsize=(11, 12))

gs = gridspec.GridSpec(
    3, 3,
    width_ratios=[1, 1, 0.05],  # última columna = colorbar
    wspace=0.05,
    hspace=0.25
)

axes = [
    [
        fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())
        for j in range(2)
    ]
    for i in range(3)
]

cax = fig.add_subplot(gs[:, 2])  # colorbar axis

# -----------------------------
# PLOT LOOP
# -----------------------------
im = None

extent = [lon[0], lon[1], lat[0], lat[1]]  # Iberia / Cantábrico

for i, row in enumerate(grid_order):
    for j, dataset_name in enumerate(row):

        ax = axes[i][j]

        dataset_resolution = datasets[dataset_name]
        path = f'{DATA_PATH_METRICS}/{dataset_name}_dana_mean_{dataset_resolution}_2024-10-29.nc'

        ax.coastlines(linewidth=0.8)
        ax.set_extent(extent, crs=ccrs.PlateCarree())

        # -------------------------
        # MISSING DATA → placeholder
        # -------------------------
        if not os.path.exists(path):
            ax.set_facecolor("white")

            ax.text(
                0.5, 0.5,
                dataset_name,
                ha='center',
                va='center',
                fontsize=11,
                transform=ax.transAxes
            )

            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(dataset_name)
            continue

        # -------------------------
        # LOAD DATA
        # -------------------------
        ds = xr.open_dataset(path)
        var_name = list(ds.data_vars)[0]
        da = ds[var_name]
        print(f"Dataset: {dataset_name}")
        print(da)

        im = ax.pcolormesh(
            ds.longitude,
            ds.latitude,
            da,
            transform=ccrs.PlateCarree(),
            shading='auto'
        )

        ax.set_title(f"{dataset_name} ({dataset_resolution})")

# -----------------------------
# COLORBAR (NO OVERLAP)
# -----------------------------
if im is not None:
    fig.colorbar(im, cax=cax)

# 👉 Guardar figura
output_path = f"{FIGURES_PATH}/comparisson_daily_mean_dana.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# 👉 Cerrar figura (importante en scripts/HPC)
plt.close(fig)

time_end = time.time()
time_elapsed = time_end - time_start
print(f"Time elapsed in minutes: {time_elapsed / 60:.2f}")
