import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#########################################################################
#REQUIRED USER INPUTS: 

#Line 19: Define file path and name for the region-filtered block ID file
#Line 20: Define file path and name for the detrended SST anomaly file
#Line 74: Define file path and name for composite output files
#Line 114: Define file path for the composited date list file
#Line 158: Define file path for the composite plot output file
#########################################################################

# Call in input files (region-filtered Block IDs and Detrended SST anomalies)
block_file = '/file_path/regionblock.yyyymm1_yyyymm2.nc'
sst_file = '/file_path/sst_anomalies_detrended.yyyymm1_yyyymm2.nc'

# Designate days to lag SST anomalies from blocked days (e.g. 5 = 5 days before blocking onset)
lag_days_list = [0, 5, 10, 15]

# Load datasets
block_ds = xr.open_dataset(block_file)
sst_ds = xr.open_dataset(sst_file)

block_obj = block_ds['object_id']
sst_anom = sst_ds['sst_anom']

# Identify blockCompDates (any nonzero point on grid)
region_nonzero = (block_obj != 0)
block_days_mask = region_nonzero.any(dim=['latitude', 'longitude'])

# Extract dates
block_times = pd.to_datetime(block_ds['time'].values)
blockCompDates = block_times[block_days_mask.values]
nonCompDates = block_times[~block_days_mask.values]

# Align to SST times
sst_times = pd.to_datetime(sst_ds['time'].values)

# Prepare storage for composites and matched dates
composites = {}
composite_dates = {}

for lag_days in lag_days_list:
    # Apply lag
    blockCompDates_lagged = blockCompDates
    nonCompDates_lagged = nonCompDates - pd.to_timedelta(lag_days, unit='D')

    # Find matching indices in SST
    blockCompIndices = sst_times.isin(blockCompDates_lagged)
    nonCompIndices = sst_times.isin(nonCompDates_lagged)

    if blockCompIndices.sum() == 0:
        print(f"No blockComp dates matched in SST dataset for lag {lag_days}!")
        continue
    if nonCompIndices.sum() == 0:
        print(f"No nonComp dates matched in SST dataset for lag {lag_days}!")
        continue

    # Extract SST composites
    blockComp = sst_anom.sel(time=sst_ds['time'][blockCompIndices]).mean(dim='time')
    nonComp = sst_anom.sel(time=sst_ds['time'][nonCompIndices]).mean(dim='time')

    # Compute difference
    compDiff = blockComp - nonComp

    # Save only compDiff to NetCDF
    out_file = f'/file_path/sst_compDiff_{lag_days}d.nc'
    compDiff.to_dataset(name='compDiff').to_netcdf(out_file)
    print(f"compDiff saved to {out_file}")

    # Adjust longitude if needed (convert 0–360 to -180–180)
    if compDiff.longitude.max() > 180:
        compDiff = compDiff.assign_coords(
            longitude=(((compDiff.longitude + 180) % 360) - 180)
        )

    # Select the desired region
    region_compDiff = compDiff.sel(
        latitude=slice(80, 35),  # descending because latitude usually goes N to S
        longitude=slice(-70, 0)  # 70°W to 0°E (negative values)
    )

    # Sort coordinates
    region_compDiff_sorted = region_compDiff.sortby(['latitude', 'longitude'])

    # Store for plotting
    composites[lag_days] = region_compDiff_sorted

    # Store matched dates
    matched_block_dates = sst_times[blockCompIndices].to_pydatetime()
    matched_nonblock_dates = sst_times[nonCompIndices].to_pydatetime()

    composite_dates[lag_days] = {
        'blockCompDates': matched_block_dates,
        'nonCompDates': matched_nonblock_dates
    }

# Save all composite dates to CSV
date_records = []
for lag, dates in composite_dates.items():
    for d in dates['blockCompDates']:
        date_records.append({'lag_days': lag, 'type': 'blockComp', 'date': d})
    for d in dates['nonCompDates']:
        date_records.append({'lag_days': lag, 'type': 'nonComp', 'date': d})

# Convert to DataFrame and save
date_df = pd.DataFrame(date_records)
date_df.sort_values(['lag_days', 'type', 'date'], inplace=True)
date_df.to_csv('/file_path/composite_dates.csv', index=False)
print("Composite date list saved to composite_dates.csv")

# Determine shared color scale with symmetric bounds around zero
all_min = min(composites[lag].min().item() for lag in composites)
all_max = max(composites[lag].max().item() for lag in composites)
absmax = max(abs(all_min), abs(all_max))
vmin, vmax = -absmax, absmax

# Create 4-panel figure
fig, axs = plt.subplots(2, 2, figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
axs = axs.flatten()

for i, lag in enumerate(lag_days_list):
    if lag not in composites:
        continue  # skip if no data for this lag
    ax = axs[i]
    im = composites[lag].plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        add_colorbar=False
    )
    ax.set_title(f"{lag}-Day Lag", fontsize=14)
    ax.coastlines(resolution='50m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.set_extent([-70, 0, 35, 80], crs=ccrs.PlateCarree())
    ax.set_xticks(np.arange(-70, 1, 10), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(35, 81, 5), crs=ccrs.PlateCarree())
    ax.set_xlabel('Longitude', fontsize=10)
    ax.set_ylabel('Latitude', fontsize=10)
    ax.tick_params(labelsize=8)
    ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', linestyle='--')

# Shared colorbar below all panels
cbar_ax = fig.add_axes([0.125, 0.1, 0.75, 0.03])  # [left, bottom, width, height]
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label('SST Composite Difference', fontsize=12)

fig.suptitle("Lagged SST Composite Differences (Blocked - Non-Blocked)", fontsize=18, fontweight='bold')
fig.subplots_adjust(wspace=0.25, hspace=0, top=0.95, bottom=0.15)

plot_file = '/file_path/sst_composite_diff.yyyymm1_yyyymm2.png'
plt.savefig(plot_file, dpi=300)
plt.close()
print(f"4-panel composite difference plot saved to {plot_file}")
