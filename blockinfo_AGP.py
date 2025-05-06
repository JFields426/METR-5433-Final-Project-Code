import numpy as np
import netCDF4
import pandas as pd
from datetime import datetime, timedelta

# Open NetCDF files
nc = netCDF4.Dataset('/share/data1/Students/jfields/finalproj/hgt.20202022_regionblock.nc', 'r')
hgt_nc = netCDF4.Dataset('/share/data1/Students/jfields/TempestExtremes/ERA5_hgt_500mb.19502022.nc', 'r')

# Read time variable and convert from "hours since 1900-01-01"
time_var = nc.variables['time']
time_origin = datetime(1900, 1, 1)
times = np.array([time_origin + timedelta(days=float(t)) for t in time_var[:]])

# Determine time resolution in days
time_resolution_days = 1

# Read lat/lon and compute area weights in steradians
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
dlat = np.radians(np.abs(lat[1] - lat[0]))
dlon = np.radians(np.abs(lon[1] - lon[0]))
lat_rad = np.radians(lat)
area_weights = np.outer(np.cos(lat_rad), np.ones(len(lon))) * dlat * dlon  # shape (lat, lon)

# Earth radius in km
R = 6371

# Open hgt data
hgt_var = hgt_nc.variables['hgt']
hgt_fill = hgt_var._FillValue

# Dictionary to store block data
block_data = {}

# Loop through time steps
for t in range(nc.dimensions['time'].size):
    block_slice = nc.variables['object_id'][t, :, :]
    valid = ~np.isnan(block_slice)
    unique_blocks = np.unique(block_slice[valid])
    
    for block in unique_blocks:
        if block == 0:
            continue  # Skip background
        
        if block not in block_data:
            block_data[block] = {
                'start_time': times[t],
                'end_time': times[t],
                'duration_steps': 0,
                'max_spatial_area_km2': 0,
                'max_size_time': None,
                'min_spatial_area_km2': np.inf,
                'min_size_time': None,
                'center_latitude': None,
                'center_longitude': None,
                'max_avg_hgt': None,
                'max_avg_hgt_time': None
            }
        
        indices = np.where(block_slice == block)
        area_km2 = np.sum(area_weights[indices]) * R**2
        
        # Update max area
        if area_km2 > block_data[block]['max_spatial_area_km2']:
            block_data[block]['max_spatial_area_km2'] = area_km2
            block_data[block]['max_size_time'] = times[t]
            block_data[block]['center_latitude'] = np.mean(lat[indices[0]])
            block_data[block]['center_longitude'] = np.mean(lon[indices[1]])
        
        # Update min area
        if area_km2 < block_data[block]['min_spatial_area_km2']:
            block_data[block]['min_spatial_area_km2'] = area_km2
            block_data[block]['min_size_time'] = times[t]
        
        # Time tracking
        block_data[block]['end_time'] = times[t]
        block_data[block]['duration_steps'] += 1
        
        # Get hgt values at the block's gridpoints
        hgt_slice = hgt_var[t, :, :]
        hgt_vals = hgt_slice[indices].astype(np.float32)
        hgt_vals[hgt_vals == hgt_fill] = np.nan
        
        avg_hgt = np.nanmean(hgt_vals)
        
        # Update max average hgt
        if block_data[block]['max_avg_hgt'] is None or avg_hgt > block_data[block]['max_avg_hgt']:
            block_data[block]['max_avg_hgt'] = avg_hgt
            block_data[block]['max_avg_hgt_time'] = times[t]

# Add duration_days field
for block in block_data:
    steps = block_data[block]['duration_steps']
    block_data[block]['duration_days'] = steps * time_resolution_days

# Convert to DataFrame with desired column order
df = pd.DataFrame.from_dict(block_data, orient='index')
df.index.name = 'block_id'
df = df.reset_index()  # make block_id a column

# Reorder and filter columns
df = df[
    ['block_id', 'start_time', 'end_time', 'duration_days', 'duration_steps',
     'min_spatial_area_km2', 'min_size_time', 'max_spatial_area_km2',
     'center_latitude', 'center_longitude',
     'max_avg_hgt', 'max_avg_hgt_time']
]

# Save CSV
output_path = '/share/data1/Students/jfields/finalproj/20202022.regionblock_summary.csv'
df.to_csv(output_path, index=False)
print(f"CSV file '{output_path}' has been created.")

# Close NetCDF files
nc.close()
hgt_nc.close()
