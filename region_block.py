import xarray as xr
import numpy as np

###################################################################
#REQUIRED USER INPUTS: 

#Line 12: Define preferred file path for block ID file
#Line 66: Define preferred file path for region-filtered block file
###################################################################

# Load the blockid data containing unique blocks
infile = '/file_path/blockid.yyyymm1_yyyymm2.nc'
ds = xr.open_dataset(infile)

object_id = ds['object_id']
lat = ds['latitude'].values
lon = ds['longitude'].values
time = ds['time'].values

# Get all unique block IDs (excluding 0, assuming 0 means no block)
unique_ids = np.unique(object_id.values)
unique_ids = unique_ids[unique_ids != 0]

# Prepare a mask to keep only selected blocks
mask = xr.zeros_like(object_id)

for block_id in unique_ids:
    # Find where this block appears
    block_mask = object_id == block_id
    # Find the first timestep where it appears
    time_sum = block_mask.sum(dim=['latitude', 'longitude'])
    first_time_idx = np.where(time_sum > 0)[0][0]

    # Get the lat/lon grid points where this block exists at first timestep
    block_points = block_mask.isel(time=first_time_idx)
    lat_points = lat[block_points.any(dim='longitude')]
    lon_points = lon[block_points.any(dim='latitude')]

    # Get indices of these points
    lat_idxs, lon_idxs = np.where(block_points.values)

    # Calculate mean lat/lon (center)
    block_lats = lat[lat_idxs]
    block_lons = lon[lon_idxs]
    center_lat = np.mean(block_lats)
    center_lon = np.mean(block_lons)
    
    # Adjust longitude to 0–360 if needed
    center_lon = center_lon % 360

    # Check if center is within 30°N–80°N, 70°W–0°E (i.e., 290°–360°)
    if 30 <= center_lat <= 80 and (290 <= center_lon <= 360 or 0 <= center_lon <= 0):
        # Keep this block
        mask = mask.where(~block_mask, block_mask)

# Set all unselected blocks to 0
filtered_object_id = object_id.where(mask, 0)

# Create a new dataset
filtered_ds = xr.Dataset(
    {'object_id': filtered_object_id},
    coords={'time': ds['time'], 'latitude': ds['latitude'], 'longitude': ds['longitude']}
)

# Save to new NetCDF file
outfile = '/file_path/regionblock.yyyymm1_yyyymm2.nc'
filtered_ds.to_netcdf(outfile)

print(f"Filtered file saved to {outfile}")
