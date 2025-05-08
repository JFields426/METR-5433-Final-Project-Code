import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

##########################################################################
#REQUIRED USER INPUTS: 

#Line 22: Define file path for the input composite files
#Line 30: Define file path and name for detrended SST anomaly input file
#Line 35: Define file path for the input composite date list file
#Line 49: [OPTIONAL] Input coordinate value for reproducablity
#Line 127: Define file path for the histogram plots output file
#Line 173: Define file path for the field-mean histogram plots output file
#Line 225: Define file path for the stippled composite plots output file
##########################################################################

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
sst_ds = xr.open_dataset('/file_path/sst_anomalies_detrended_standardized.yyyymm1_yyyymm2.nc')
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
sst_mean = sst_mean.sel(latitude=slice(80, 35), longitude=slice(-70, 0)).sortby(['latitude', 'longitude'])
ocean_mask = ~np.isnan(sst_mean)

# Manual override for reproducibility
manual_lat_lon = None  # Example: (50.0, -40.0)

if manual_lat_lon:
    lat_val, lon_val = manual_lat_lon
else:
    valid_lat, valid_lon = np.where(ocean_mask.values)
    chosen_idx = np.random.choice(len(valid_lat))
    lat_val = sst_mean.latitude.values[valid_lat[chosen_idx]]
    lon_val = sst_mean.longitude.values[valid_lon[chosen_idx]]

print(f"Using ocean point at lat={lat_val:.2f}, lon={lon_val:.2f}")

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
    pvals = np.mean(np.abs(boot_array) >= np.abs(actual_diff), axis=0)
    pval_da = xr.DataArray(pvals, coords=composites[lag].coords, dims=composites[lag].dims)
    bootstrap_pvals[lag] = pval_da

# Plot histograms of sampled differences at selected point
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

for i, lag in enumerate(lags):
    boot_array = boot_diffs_all[lag]

    # Extract values at selected lat/lon using xarray indexing
    sample_vals = []
    for arr in boot_array:
        da = xr.DataArray(arr, coords=composites[lag].coords, dims=composites[lag].dims)
        val = da.sel(latitude=lat_val, longitude=lon_val, method="nearest").item()
        sample_vals.append(val)
    sample_vals = np.array(sample_vals)

    actual_val = composites[lag].sel(latitude=lat_val, longitude=lon_val, method="nearest").item()

    p5, p95 = np.percentile(sample_vals, [2.5, 97.5])
    p10, p90 = np.percentile(sample_vals, [5, 95])
    
    x_pad = (p95 - p5) * 0.3  # extra space on either side
    x_min = min(p5, actual_val) - x_pad
    x_max = max(p95, actual_val) + x_pad

    ax = axs[i]
    ax.hist(sample_vals, bins=60, density=True, color='skyblue', edgecolor='gray')
    ax.axvline(p5, color='red', linestyle='--', label='p<0.05')
    ax.axvline(p95, color='red', linestyle='--')
    ax.axvline(p10, color='orange', linestyle='--', label='p<0.10')
    ax.axvline(p90, color='orange', linestyle='--')
    ax.axvline(actual_val, color='black', linestyle='-', label='Actual Value')
    ax.set_xlim(x_min, x_max)
    ax.set_title(f'{lag}-Day Lag at ({lat_val:.1f}N, {lon_val:.1f}E)', fontsize=12)
    ax.set_xlabel('Composite Difference (σ)')
    ax.set_ylabel('PDF')
    ax.legend(fontsize=8)
    ax.grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.savefig('/file_path/bootstrap_histograms.yyyymm1_yyyym2.png', dpi=300)
plt.close()

# Plot histograms of sampled differences for the field mean

fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

for i, lag in enumerate(lags):
    boot_array = boot_diffs_all[lag]
    # Collapse each bootstrap sample into a spatial mean over the full field
    field_means = []
    for arr in boot_array:
        da = xr.DataArray(arr, coords=composites[lag].coords, dims=composites[lag].dims)
        # Apply ocean mask to only keep valid ocean points
        masked_da = da.where(ocean_mask)
        mean_val = masked_da.mean().item()
        field_means.append(mean_val)
    field_means = np.array(field_means)

    # Actual composite mean
    actual_field_mean = composites[lag].where(ocean_mask).mean().item()

    # Percentiles
    p5, p95 = np.percentile(field_means, [2.5, 97.5])
    p10, p90 = np.percentile(field_means, [5, 95])

    x_pad = (p95 - p5) * 0.3
    x_min = min(p5, actual_field_mean) - x_pad
    x_max = max(p95, actual_field_mean) + x_pad

    ax = axs[i]
    ax.hist(field_means, bins=60, density=True, color='lightgreen', edgecolor='gray')
    ax.axvline(p5, color='red', linestyle='--', label='p<0.05')
    ax.axvline(p95, color='red', linestyle='--')
    ax.axvline(p10, color='orange', linestyle='--', label='p<0.10')
    ax.axvline(p90, color='orange', linestyle='--')
    ax.axvline(actual_field_mean, color='black', linestyle='-', label='Actual Field Mean')
    ax.set_xlim(x_min, x_max)
    ax.set_title(f'{lag}-Day Lag Field Mean', fontsize=12)
    ax.set_xlabel('Field Composite Difference Mean (σ)')
    ax.set_ylabel('PDF')
    ax.legend(fontsize=8)
    ax.grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.savefig('/file_path/bootstrap_field_histograms.yyyymm1_yyyymm2png', dpi=300)
plt.close()

# Plot Composite Maps
all_min = min(composites[lag].min().item() for lag in composites)
all_max = max(composites[lag].max().item() for lag in composites)
absmax = max(abs(all_min), abs(all_max))
vmin, vmax = -absmax, absmax

fig, axs = plt.subplots(2, 2, figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
axs = axs.flatten()

for i, lag in enumerate(lags):
    comp = composites[lag]
    pval = bootstrap_pvals[lag]

    sst_mean_comp = sst_anom.mean(dim='time')
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

    ax.set_title(f"{lag}-Day Lag", fontsize=12)
    ax.coastlines(resolution='50m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.set_extent([-70, 0, 35, 80], crs=ccrs.PlateCarree())
    ax.set_xticks(np.arange(-70, 1, 10), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(35, 81, 5), crs=ccrs.PlateCarree())
    ax.set_xlabel('Longitude', fontsize=10)
    ax.set_ylabel('Latitude', fontsize=10)
    ax.tick_params(labelsize=8)
    ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', linestyle='--')

cbar_ax = fig.add_axes([0.125, 0.1, 0.75, 0.03])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label('SST Composite Difference (σ)')

fig.suptitle("Lagged Standardized SST Anomaly Composite Differences\n(Blocked - Non-Blocked) (p < 0.05)", fontsize=18, fontweight='bold')
fig.subplots_adjust(wspace=0.25, hspace=0, top=0.95, bottom=0.15)

plt.savefig('/file_path/sst_composite.yyyymm1_yyyymm2.png', dpi=300)
plt.close()
