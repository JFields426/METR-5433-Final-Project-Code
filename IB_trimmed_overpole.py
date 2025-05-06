import xarray as xr
import numpy as np
import glob

# Define the latitude ranges, including adjustments for latitudes above 70 degrees
lat0_range = np.arange(35, 90.5, 0.5)  # Full range from 35 to 90 degrees
latS_range = np.where(lat0_range > 70, lat0_range - (90 - lat0_range), lat0_range - 20)
latN_range = np.where(lat0_range > 70, lat0_range + (90 - lat0_range), lat0_range + 20)
lat15S_range = np.where(lat0_range > 70, lat0_range - (90 - lat0_range) * 0.75, lat0_range - 15)
lat30S_range = np.where(lat0_range > 70, lat0_range - (90 - lat0_range) * 1.5, lat0_range - 30)

pressure_level = 500

# Define the file path pattern
file_path_pattern = '/data/deluge/reanalysis/REANALYSIS/ERA5/3D/4xdaily/hgt/hgt.*.nc'

# Get a list of all files matching the pattern
all_files = glob.glob(file_path_pattern)

# Filter files based on the specified range (1950-01 to 2022-12)
filtered_files = [f for f in all_files if 'hgt.202001' <= f.split('/')[-1] <= 'hgt.202212.nc']

# Ensure filtered_files is sorted correctly
filtered_files.sort(key=lambda x: int(x.split('/')[-1].split('.')[1]))

# Initialize lists to store results
valid_masks = []
GHGS_list = []
GHGN_list = []
times = []

# Process each file
for file_path in filtered_files:

    # Load the file with automatic chunking
    ds = xr.open_dataset(file_path, chunks='auto')

    # Select data for the specific pressure level
    ds_500 = ds.sel(level=pressure_level)

    # Resample to daily means
    ds_daily = ds_500.resample(time='1D').mean().persist()

    # Drop 'level' since it's now constant
    ds_daily = ds_daily.drop_vars('level')

    # Select latitude ranges
    Z500_lat0 = ds_daily['hgt'].sel(latitude=lat0_range, method='nearest')
    Z500_latS = ds_daily['hgt'].sel(latitude=latS_range, method='nearest')
    Z500_latN = ds_daily['hgt'].sel(latitude=latN_range, method='nearest')
    Z500_lat15S = ds_daily['hgt'].sel(latitude=lat15S_range, method='nearest')
    Z500_lat30S = ds_daily['hgt'].sel(latitude=lat30S_range, method='nearest')

    # Calculate GHGS and GHGN
    GHGS = (Z500_lat0.values - Z500_latS.values) / 20
    GHGN = (Z500_latN.values - Z500_lat0.values) / 20

    # Calculate the additional gradient condition
    gradient_S = (Z500_lat15S.values - Z500_lat30S.values)

    # Create a mask where all conditions are met
    valid_mask = ((GHGS > 0) & (GHGN < -10) & (gradient_S < 0)).astype(int)

    # Create valid_mask as a DataArray with the same coordinates as Z500_lat0
    valid_mask = xr.DataArray(valid_mask, coords=Z500_lat0.coords, dims=Z500_lat0.dims)

    # Append results to lists
    valid_masks.append(valid_mask)
    GHGS_list.append(xr.DataArray(GHGS, coords=Z500_lat0.coords, dims=Z500_lat0.dims))
    GHGN_list.append(xr.DataArray(GHGN, coords=Z500_lat0.coords, dims=Z500_lat0.dims))
    times.append(Z500_lat0['time'])
    print(file_path)

# Combine all the results into a single Dataset
valid_mask_combined = xr.concat(valid_masks, dim='time')
GHGS_combined = xr.concat(GHGS_list, dim='time')
GHGN_combined = xr.concat(GHGN_list, dim='time')

# Apply a rolling window of 5 days to enforce 5-day persistence
rolling_block = valid_mask_combined.rolling(time=5, center=True).sum()
persistent_block = (rolling_block >= 5).astype(int)
valid_mask_combined = persistent_block

# Create the final result dataset
result = xr.Dataset({
    'block_tag': (('time', 'latitude', 'longitude'), valid_mask_combined.data),
    'GHGS': (('time', 'latitude', 'longitude'), GHGS_combined.data),
    'GHGN': (('time', 'latitude', 'longitude'), GHGN_combined.data)
}, coords={
    'time': np.concatenate(times),  # Combine all time coordinates
    'latitude': Z500_lat0['latitude'],
    'longitude': Z500_lat0['longitude']
})

# Set encodings for the nc file
encoding = {
    'block_tag': {"dtype": "float64", "complevel": 9, "zlib": True},
    'GHGS': {"dtype": "float64", "complevel": 9, "zlib": True},
    'GHGN': {"dtype": "float64", "complevel": 9, "zlib": True},
    'time': {"dtype": "int32", "complevel": 9, "zlib": True, "calendar": "gregorian", "units": "days since 1900-01-01"},
    'latitude': {"dtype": "float64", "complevel": 9, "zlib": True},
    'longitude': {"dtype": "float64", "complevel": 9, "zlib": True},
}

# Save the dataset as a NetCDF file
output_path = '/share/data1/Students/jfields/BlockingDataset/block_tag.20202022.nc'
result.to_netcdf(output_path, format='NETCDF4_CLASSIC', encoding=encoding)
