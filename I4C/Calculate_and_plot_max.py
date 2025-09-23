import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import glob
import os
import sys

matplotlib.use("Agg")

# -------------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------------
if len(sys.argv) < 2:
    print("Usage: python code.py <freq> [varname]")
    sys.exit(1)

freq = sys.argv[1]  # required
varname = sys.argv[2] if len(sys.argv) > 2 else "pr"  # optional, defaults to "pr"

print("Frequency:", freq)
print("Variable name:", varname)

# -------------------------------------------------------------------
# Institutions and domains
# -------------------------------------------------------------------
institutions = ["BCCR-UCAN", "BCCR-UCAN", "CNRM-MF", "ICTP"]
domains = ["EUR-12", "ALPX-3", "ALPX-3", "VAL-3"]

# -------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------
def calculate_yearly_max(
    path_data,
    period,
    varname="pr",
    lon_city=-0.375, lat_city=39.467,
    n_grid=30,
    freq="1hr"
):
    """
    Calculate mean yearly maximum and absolute maximum precipitation
    around a city using high-resolution climate data.
    """

    # Load files
    all_files = sorted(glob.glob(f"{path_data}/{freq}/{varname}/*/*.nc"))
    if period == "historical":
        files = [f for f in all_files if any(f"_{year}" in f for year in range(1995, 2015))]
    elif period == "ssp370":
        files = all_files[-20:]
    else:
        files = all_files  # fallback: keep everything

    if not files:
        raise FileNotFoundError(f"No files found for {path_data} / {period}")
    else:
        print(f"{len(files)} files found for {period} at {path_data}")

    ds = xr.open_mfdataset(files, combine="nested", concat_dim="time")

    # Convert units if needed
    if ds[varname].attrs.get("units", "") == "kg m-2 s-1":
        if freq == "1hr":
            ds[varname] *= 3600
            ds[varname].attrs["units"] = "mm/1hr"
        elif freq in ["day", "daily"]:
            ds[varname] *= 86400
            ds[varname].attrs["units"] = "mm/day"

    # Find closest grid point to the city
    lat = ds["lat"].values
    lon = ds["lon"].values
    abs_diff = np.abs(lat - lat_city) + np.abs(lon - lon_city)
    iy, ix = np.unravel_index(np.argmin(abs_diff), lat.shape)

    # Slice area around the city
    half_box = n_grid // 2
    y_slice = slice(max(0, iy - half_box), iy + half_box)
    x_slice = slice(max(0, ix - half_box), ix + half_box)

    # Subset and compute yearly max
    ds_subset = ds[varname].isel(y=y_slice, x=x_slice)
    ds_yearly_max = ds_subset.groupby("time.year").max("time")
    ds_mean_max = ds_yearly_max.mean(dim="year").compute()

    # Compute absolute maximum over the whole period
    ds_max = ds_subset.max("time").compute()

    # Load land-sea mask
    mask_files = sorted(glob.glob(f"{path_data}/fx/sftlf/*/sftlf_*.nc"))
    mask = xr.open_mfdataset(mask_files)["sftlf"].isel(y=y_slice, x=x_slice)

    return ds_mean_max, ds_max, mask


def plot_data(
    datain,
    mask,
    domain,
    institution,
    period,
    varname="pr",
    lon_city=-0.375, lat_city=39.467,
    n_grid=30,
    vmin=0,
    vmax=70,
    cmap="viridis",
    name="max"
):
    """Plot a single map of data with mask and city marker."""

    fig, ax = plt.subplots(figsize=(8, 6))

    # Main data plot
    datain.plot(
        ax=ax,
        x="lon", y="lat",
        cmap=cmap,
        vmin=vmin, vmax=vmax,
        cbar_kwargs={"label": "Precipitation (mm)"}
    )

    # Contour for land-sea mask
    mask.plot.contour(
        ax=ax,
        x="lon", y="lat",
        levels=[50],  # 50% land threshold
        colors="black",
        linewidths=1.0,
        zorder=10
    )

    # City marker
    ax.plot(lon_city, lat_city,
            marker="o", markersize=10,
            markerfacecolor="none", markeredgecolor="red", zorder=11)

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"{name} {varname} ({institution}, {domain}, {period})", fontsize=12)
    plt.tight_layout()

    # Save figure
    outdir = f"Figures/{institution}"
    os.makedirs(outdir, exist_ok=True)
    outfile = f"{outdir}/{name}_{varname}_{domain}_{institution}_{period}.png"
    plt.savefig(outfile, dpi=300)
    plt.close(fig)  # don't fill memory
    print(f"Figure saved to: {outfile}")


# -------------------------------------------------------------------
# Custom colormap
# -------------------------------------------------------------------
colors = [
    "#ffffff", "#e6f7ff", "#cceeff", "#b3e6ff", "#99ddff", "#80d4ff",
    "#66ccff", "#4dc3ff", "#33bbff", "#1ab2ff", "#00aaff",
    "#00bfb3", "#00cc88", "#66d96d", "#c2e255", "#ffe945",
    "#ffd032", "#ffac24", "#ff7b22", "#d73027"
]
bounds = np.arange(0, 21, 1)  # 21 boundaries → 20 bins
cmap = ListedColormap(colors)
norm = BoundaryNorm(bounds, ncolors=len(colors))

# -------------------------------------------------------------------
# Historical runs
# -------------------------------------------------------------------
period = "historical"
results_histo = {}

for institution, domain in zip(institutions, domains):
    print(institution, domain)
    root = "/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/"
    if institution == "BCCR-UCAN":
        path_data = f"{root}/{domain}/{institution}/NorESM2-MM/{period}/r1i1p1f1/WRF451I/v1-r1/"
    elif institution == "CNRM-MF":
        path_data = f"{root}/{domain}/{institution}/CNRM-ESM2-1/{period}/r1i1p1f2/CNRM-AROME46t1/v1-r1/"
    elif institution == "ICTP":
        path_data = f"{root}/{domain}/{institution}/EC-Earth3-Veg/{period}/r1i1p1f1/RegCM5-0/v1-r1/"
    else:
        continue

    n_grid = 30 if domain == "EUR-12" else 90

    try:
        results_histo[(institution, domain)] = calculate_yearly_max(
            path_data, period=period, varname=varname,
            lon_city=-0.375, lat_city=39.467, n_grid=n_grid, freq=freq
        )
    except Exception as e:
        print(f"Skipping {institution} - {domain} due to error: {e}")
        continue

# Plot historical
for (institution, domain), value in results_histo.items():
    plot_data(value[0], value[2], domain, institution, period, varname, vmin=0, vmax=50, cmap=cmap, name=f"mean_max_{freq}")
    plot_data(value[1], value[2], domain, institution, period, varname, vmin=0, vmax=180, cmap=cmap, name=f"max_{freq}")

# -------------------------------------------------------------------
# SSP370 runs
# -------------------------------------------------------------------
period = "ssp370"
results_ssp = {}

for institution, domain in zip(institutions, domains):
    print(institution, domain)
    root = "/gpfs/projects/meteo/WORK/ASNA/projects/impetus/DATA/I4C/CMIP6/DD/"
    if institution == "BCCR-UCAN":
        path_data = f"{root}/{domain}/{institution}/NorESM2-MM/{period}/r1i1p1f1/WRF451I/v1-r1/"
    elif institution == "CNRM-MF":
        path_data = f"{root}/{domain}/{institution}/CNRM-ESM2-1/{period}/r1i1p1f2/CNRM-AROME46t1/v1-r1/"
    elif institution == "ICTP":
        path_data = f"{root}/{domain}/{institution}/EC-Earth3-Veg/{period}/r1i1p1f1/RegCM5-0/v1-r1/"
    else:
        continue

    n_grid = 30 if domain == "EUR-12" else 90

    try:
        results_ssp[(institution, domain)] = calculate_yearly_max(
            path_data, period=period, varname=varname,
            lon_city=-0.375, lat_city=39.467, n_grid=n_grid, freq=freq
        )
    except Exception as e:
        print(f"Skipping {institution} - {domain} due to error: {e}")
        continue

# Plot SSP370
for (institution, domain), value in results_ssp.items():
    plot_data(value[0], value[2], domain, institution, period, varname, vmin=0, vmax=50, cmap=cmap, name=f"mean_max_{freq}")
    plot_data(value[1], value[2], domain, institution, period, varname, vmin=0, vmax=180, cmap=cmap, name=f"max_{freq}")

# -------------------------------------------------------------------
# Differences (ssp370 - historical)
# -------------------------------------------------------------------
delta = {}
for key, ssp_value in results_ssp.items():
    if key in results_histo:
        hist_value = results_histo[key]
        delta[key] = (
            ssp_value[0] - hist_value[0],  # mean_max diff
            ssp_value[1] - hist_value[1],  # max diff
            ssp_value[2]                   # keep mask from ssp
        )
    else:
        print(f"Skipping {key} as it is not in historical results")

# Plot deltas
for (institution, domain), value in delta.items():
    plot_data(value[0], results_ssp[institution, domain][2], domain, institution, period, varname, vmin=-20, vmax=20, cmap="BrBG", name=f"delta_mean_max_{freq}")
    plot_data(value[1], results_ssp[institution, domain][2], domain, institution, period, varname, vmin=-20, vmax=20, cmap="BrBG", name=f"delta_max_{freq}")
