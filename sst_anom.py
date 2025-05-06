import xarray as xr
import numpy as np
import glob
import re

#########################################################################
#REQUIRED USER INPUTS: 

#Line 16: Define file path for ERA5 SST input files
#Line 17: Define file path for SST climatology file (update file name if different time range was used)
#Line 27: Change yyyymm1 and yyyymm2 to the start and end files for your selected time range
#Line 62: Define file path and file name for the output file
#########################################################################

#Call in ERA5 sea surface temperature reanalysis data (assumes data is already downloaded)
file_path_pattern = '/file_path/sst.*.nc'
climo_path = '/file_path/sst_climo_1991_2020.nc'

#Calls in climatology file to calculate anomalies
climatology = xr.open_dataset(climo_path)['sst']

#Isolate files for the date range you want to calculate anomalies for
all_files = glob.glob(file_path_pattern)
target_files = [
    f for f in all_files
    if re.match(r'.*/sst\.(\d{6})\.nc$', f)
    and yyyymm1 <= int(re.search(r'sst\.(\d{6})\.nc$', f).group(1)) <= yyyymm2
]

target_files.sort(key=lambda x: int(re.search(r'sst\.(\d{6})\.nc$', x).group(1)))

#Define spatial domain (30N-80N, 70W-0E)
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

#Combine anomalies into one dataset and save
sst_anomalies_combined = xr.concat(anomaly_list, dim='time')
anom_output_path = '/file_path/sst_anomalies.yyyymm1_yyyymm2.nc'
sst_anomalies_combined.to_netcdf(anom_output_path, format='NETCDF4_CLASSIC')

print(f"SST anomalies saved to {anom_output_path}")
