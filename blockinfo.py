import numpy as np
import netCDF4
import pandas as pd
from datetime import datetime, timedelta

##################################################################################################################################
#REQUIRED USER INPUTS: 

#Line 15: Define file path for the region-filtered block ID file
#Line 23: [OPTIONAL] Change the time resolution to reflect that of the input file. Only do this if the time resolution is not daily.
#Line 99: Define the file path for the region-filtered block ID summary file
##################################################################################################################################

# Open NetCDF files
nc = netCDF4.Dataset('/file_path/regionblock.yyyymm1_yyyymm2.nc', 'r')

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
     'center_latitude', 'center_longitude']
]

# Save CSV
output_path = '/file_path/regionblock_summary.yyyymm1_yyyymm2.csv'
df.to_csv(output_path, index=False)
print(f"CSV file '{output_path}' has been created.")

# Close NetCDF files
nc.close()
hgt_nc.close()
