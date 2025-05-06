import xarray as xr
import numpy as np
import glob
import re

# ---------------------
# File paths and filtering
# ---------------------
file_path_pattern = '/data/deluge/reanalysis/REANALYSIS/ERA5/2D/4xdaily/sst/sst.*.nc'
all_files = glob.glob(file_path_pattern)

# Keep only files matching sst.YYYYMM.nc
filtered_files = [f for f in all_files if re.match(r'.*/sst\.\d{6}\.nc$', f)]

# Sort by extracted YYYYMM
filtered_files.sort(key=lambda x: int(re.search(r'sst\.(\d{6})\.nc$', x).group(1)))

# Further filter for climatology period
climo_files = [f for f in filtered_files if 199101 <= int(re.search(r'sst\.(\d{6})\.nc$', f).group(1)) <= 202012]

# ---------------------
# Domain selection
# ---------------------
lat_min, lat_max = 30, 80
lon_min, lon_max = (-70 % 360), (0 % 360)  # 290 to 0 (wrap-around handled)

# ---------------------
# Load and process
# ---------------------
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

    # Mask NaNs (land) for consistency
    sst_daily = sst_daily.where(np.isfinite(sst_daily))

    daily_sst_list.append(sst_daily)
    print(f"Loaded: {file_path}")

# ---------------------
# Combine and compute climatology
# ---------------------
combined_sst = xr.concat(daily_sst_list, dim='time')

# Group by day of year
climatology = combined_sst.groupby('time.dayofyear').mean('time', skipna=True)

# ---------------------
# Save
# ---------------------
climo_output_path = '/share/data1/Students/jfields/finalproj/sst_climo_1991_2020.nc'
climatology.to_netcdf(climo_output_path, format='NETCDF4_CLASSIC')
print(f"Climatology saved to {climo_output_path}")
