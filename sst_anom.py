import xarray as xr
import numpy as np
import glob
import re

# ---------------------
# Paths
# ---------------------
file_path_pattern = '/data/deluge/reanalysis/REANALYSIS/ERA5/2D/4xdaily/sst/sst.*.nc'
climo_path = '/share/data1/Students/jfields/finalproj/sst_climo_1991_2020.nc'

# ---------------------
# Load climatology
# ---------------------
climatology = xr.open_dataset(climo_path)['sst']

# ---------------------
# Load SST files 2000–2022
# ---------------------
all_files = glob.glob(file_path_pattern)
target_files = [
    f for f in all_files
    if re.match(r'.*/sst\.(\d{6})\.nc$', f)
    and 195001 <= int(re.search(r'sst\.(\d{6})\.nc$', f).group(1)) <= 202212
]

target_files.sort(key=lambda x: int(re.search(r'sst\.(\d{6})\.nc$', x).group(1)))

# ---------------------
# Domain
# ---------------------
lat_min, lat_max = 30, 80
lon_min, lon_max = (-70 % 360), (0 % 360)

anomaly_list = []

for file_path in target_files:
    ds = xr.open_dataset(file_path)
    sst = ds['sst'].sel(latitude=slice(lat_max, lat_min))

    # Handle wrap-around for longitudes
    if lon_min > lon_max:
        sst = xr.concat([
            sst.sel(longitude=slice(lon_min, 360)),
            sst.sel(longitude=slice(0, lon_max))
        ], dim='longitude')
    else:
        sst = sst.sel(longitude=slice(lon_min, lon_max))

    sst_daily = sst.resample(time='1D').mean()
    sst_daily = sst_daily.where(np.isfinite(sst_daily))

    # Subtract climatology: match by day of year
    anomalies = sst_daily.groupby('time.dayofyear') - climatology
    anomaly_list.append(anomalies)

    print(f"Processed anomalies: {file_path}")

# ---------------------
# Combine and Save
# ---------------------
sst_anomalies_combined = xr.concat(anomaly_list, dim='time')
anom_output_path = '/share/data1/Students/jfields/finalproj/sst_anomalies_1950_2022.nc'
sst_anomalies_combined.to_netcdf(anom_output_path, format='NETCDF4_CLASSIC')

print(f"SST anomalies saved to {anom_output_path}")
