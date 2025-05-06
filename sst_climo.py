import xarray as xr
import numpy as np
import glob
import re

###################################################################################################
#REQUIRED USER INPUTS: 

#Line 16: Define file path for ERA5 SST input files
#Line 21: [OPTIONAL] Change file range if using climatology other than 1991-2020
#Line 56: Define file path for output file (OPTIONAL: Change file name to reflect climatology range)
####################################################################################################

#Call in ERA5 sea surface temperature reanalysis data (assumes data is already downloaded)
file_path_pattern = '/file_path/sst.*.nc'
all_files = glob.glob(file_path_pattern)

# Isolate files for climatology time range (1991-2020 used in this case)
filtered_files = [f for f in all_files if re.match(r'.*/sst\.\d{6}\.nc$', f)]
filtered_files.sort(key=lambda x: int(re.search(r'sst\.(\d{6})\.nc$', x).group(1)))
climo_files = [f for f in filtered_files if 199101 <= int(re.search(r'sst\.(\d{6})\.nc$', f).group(1)) <= 202012]

lat_min, lat_max = 30, 80
lon_min, lon_max = (-70 % 360), (0 % 360)  # 290 to 0 (wrap-around handled)

daily_sst_list = []

for file_path in climo_files:
    ds = xr.open_dataset(file_path)
    sst = ds['sst'].sel(latitude=slice(lat_max, lat_min))

    # Handle 0–360 wrap
    if lon_min > lon_max:
        sst = xr.concat([
            sst.sel(longitude=slice(lon_min, 360)),
            sst.sel(longitude=slice(0, lon_max))
        ], dim='longitude')
    else:
        sst = sst.sel(longitude=slice(lon_min, lon_max))

    # Resample to daily mean
    sst_daily = sst.resample(time='1D').mean()

    # Mask NaNs (land) for future statistical functions
    sst_daily = sst_daily.where(np.isfinite(sst_daily))

    daily_sst_list.append(sst_daily)
    print(f"Loaded: {file_path}")

#Combine into a single dataset
combined_sst = xr.concat(daily_sst_list, dim='time')

# Group by day of year
climatology = combined_sst.groupby('time.dayofyear').mean('time', skipna=True)

climo_output_path = '/file_path/sst_climo_1991_2020.nc'
climatology.to_netcdf(climo_output_path, format='NETCDF4_CLASSIC')
print(f"Climatology saved to {climo_output_path}")
