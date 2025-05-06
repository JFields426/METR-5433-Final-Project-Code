import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#########################################################################
#REQUIRED USER INPUTS: 

#Line 22: Define file path for the input composite files
#Line 30: Define file path and name for detrended SST anomaly input file
#Line 35: Define file path for the input composite date list file
#Line 117: Define file path for the histogram plots output file
#Line 170: Define file path for the stippled composite plots output file
#########################################################################

# Load composite differences
composites = {}
lags = [0, 5, 10, 15]
for lag in lags:
    ds = xr.open_dataset(f'/file_path/sst_compDiff_{lag}d.nc')
    comp = ds['compDiff']
    if comp.longitude.max() > 180:
        comp = comp.assign_coords(longitude=(((comp.longitude + 180) % 360) - 180))
    comp = comp.sel(latitude=slice(80, 35), longitude=slice(-70, 0)).sortby(['latitude', 'longitude'])
    composites[lag] = comp

# Load SST anomalies
sst_ds = xr.open_dataset('/file_path/sst_anomalies_detrended.yyyymm1_yyyymm2.nc')
sst_anom = sst_ds['sst_anom']
sst_times = pd.to_datetime(sst_ds['time'].values)

# Load composite dates
date_df = pd.read_csv('/file_path/composite_dates.csv')
date_df['date'] = pd.to_datetime(date_df['date'])

# Prepare bootstrap results
bootstrap_pvals = {}
boot_diffs_all = {}
n_bootstrap = 5000

# Create ocean mask & select one random ocean point
sst_mean = sst_anom.mean(dim='time')
if sst_mean.longitude.max() > 180:
    sst_mean = sst_mean.assign_coords(longitude=(((sst_mean.longitude + 180) % 360) - 180))
sst_mean = sst_mean.sel(latitude=slice(80, 35), longitude=slice(-70, 0)).sortby(['latitude', 'longitude'])

ocean_mask = ~np.isnan(sst_mean)
valid_points = np.argwhere(ocean_mask.values)

# Optional manual override (row, col indices)
manual_row_col = None 
if manual_row_col:
    random_row, random_col = manual_row_col
else:
    chosen_idx = np.random.choice(len(valid_points))
    random_row, random_col = tuple(valid_points[chosen_idx])

lat_val = sst_mean.latitude.values[random_row]
lon_val = sst_mean.longitude.values[random_col]
print(f"Using ocean point at (row={random_row}, col={random_col}) -> lat={lat_val:.2f}, lon={lon_val:.2f}")

# Iterate through each lag
for lag in lags:
    lag_df = date_df[date_df['lag_days'] == lag]
    block_dates = lag_df[lag_df['type'] == 'blockComp']['date'].values
    nonblock_dates = lag_df[lag_df['type'] == 'nonComp']['date'].values
    n_block = len(block_dates)
    n_nonblock = len(nonblock_dates)

    block_idx = np.where(sst_times.isin(block_dates))[0]
    nonblock_idx = np.where(sst_times.isin(nonblock_dates))[0]

    # Bootstrap differences
    boot_diffs = []
    for _ in range(n_bootstrap):
        block_sample = sst_anom.isel(time=np.random.choice(block_idx, size=n_block, replace=True)).mean(dim='time')
        nonblock_sample = sst_anom.isel(time=np.random.choice(nonblock_idx, size=n_nonblock, replace=True)).mean(dim='time')
        boot_diffs.append((block_sample - nonblock_sample).values)
    boot_array = np.stack(boot_diffs, axis=0)
    boot_diffs_all[lag] = boot_array

    # Calculate empirical two-tailed p-values
    actual_diff = composites[lag].values
    pvals = np.mean(np.abs(boot_array - np.mean(boot_array, axis=0)) >= np.abs(actual_diff), axis=0)
    pval_da = xr.DataArray(pvals, coords=composites[lag].coords, dims=composites[lag].dims)
    bootstrap_pvals[lag] = pval_da

# Plot histograms of sampled differences at selected point
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

for i, lag in enumerate(lags):
    boot_array = boot_diffs_all[lag]
    sample_vals = [boot[random_row, random_col] for boot in boot_array]
    actual_val = composites[lag].values[random_row, random_col]

    p5, p95 = np.percentile(sample_vals, [2.5, 97.5])
    p10, p90 = np.percentile(sample_vals, [5, 95])
    x_max = max(abs(min(sample_vals)), abs(max(sample_vals)), abs(actual_val)) * 1.1

    ax = axs[i]
    ax.hist(sample_vals, bins=50, density=True, color='skyblue', edgecolor='gray')
    ax.axvline(p5, color='red', linestyle='--', label='p<0.05')
    ax.axvline(p95, color='red', linestyle='--')
    ax.axvline(p10, color='orange', linestyle='--', label='p<0.10')
    ax.axvline(p90, color='orange', linestyle='--')
    ax.axvline(actual_val, color='black', linestyle='-', label='Actual Value')
    ax.set_xlim(-x_max, x_max)
    ax.set_title(f'Lag {lag}d @ ({lat_val:.1f}N, {lon_val:.1f}E)')
    ax.set_xlabel('Composite Difference')
    ax.set_ylabel('PDF')
    ax.legend()

plt.tight_layout()
plt.savefig('/file_path/bootstrap_histograms_random_point.png', dpi=300)
plt.close()

# Plot composite maps with stippling only over ocean
vmin = min([composites[lag].min().item() for lag in composites])
vmax = max([composites[lag].max().item() for lag in composites])

fig, axs = plt.subplots(2, 2, figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
axs = axs.flatten()

for i, lag in enumerate(lags):
    comp = composites[lag]
    pval = bootstrap_pvals[lag]

    # Recompute ocean mask for this composite
    sst_mean_comp = sst_anom.mean(dim='time')
    if sst_mean_comp.longitude.max() > 180:
        sst_mean_comp = sst_mean_comp.assign_coords(longitude=(((sst_mean_comp.longitude + 180) % 360) - 180))
    sst_mean_comp = sst_mean_comp.sel(latitude=comp.latitude, longitude=comp.longitude).sortby(['latitude', 'longitude'])
    ocean_mask = ~np.isnan(sst_mean_comp)

    sig_mask = (pval < 0.05) & ocean_mask

    ax = axs[i]
    im = comp.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        add_colorbar=False
    )
    lons, lats = np.meshgrid(comp.longitude, comp.latitude)
    ax.plot(lons[sig_mask], lats[sig_mask], 'k.', markersize=0.5, transform=ccrs.PlateCarree())

    ax.set_title(f"Lag {lag} days", fontsize=12)
    ax.coastlines(resolution='50m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.set_extent([-70, 0, 35, 80], crs=ccrs.PlateCarree())
    ax.set_xticks(np.arange(-70, 1, 10), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(35, 81, 10), crs=ccrs.PlateCarree())
    ax.set_xlabel('Longitude', fontsize=10)
    ax.set_ylabel('Latitude', fontsize=10)
    ax.tick_params(labelsize=8)
    ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', linestyle='--')

cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label('SST Composite Difference')

fig.suptitle("SST Composite Difference with Significant Ocean Points (p < 0.05)\n35°N–80°N, 70°W–0°E", fontsize=16)
fig.subplots_adjust(wspace=0.1, hspace=0.2, top=0.92, bottom=0.15)

plt.savefig('/file_path/sst_composite_diff_with_significance.png', dpi=300)
plt.close()
