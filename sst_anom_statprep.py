import xarray as xr

# ---------------------
# Load anomalies
# ---------------------
anom_input_path = '/share/data1/Students/jfields/finalproj/sst_anomalies_2020_2022.nc'
anomalies = xr.open_dataset(anom_input_path)['sst']

# ---------------------
# Detrend anomalies (remove linear trend at each grid point)
# ---------------------
# Add numeric time coordinate
time_numeric = xr.DataArray(
    anomalies['time'].dt.year + (anomalies['time'].dt.dayofyear - 1) / 365.25,
    dims='time'
)

# Fit a 1st degree polynomial (linear) across time
coeffs = anomalies.polyfit(dim='time', deg=1, skipna=True)

# Evaluate the trend at each time step
trend = xr.polyval(time_numeric, coeffs.polyfit_coefficients)

# Remove trend from original anomalies
detrended_anomalies = anomalies - trend

detrended_anomalies.name = 'sst_anom'

# ---------------------
# Standardize anomalies (mean 0, std 1 at each grid point)
# ---------------------
#mean = detrended_anomalies.mean(dim='time', skipna=True)
#std = detrended_anomalies.std(dim='time', skipna=True)

#standardized_anomalies = (detrended_anomalies - mean) / std

# Name the variable explicitly
#standardized_anomalies.name = 'sst_anom_std'

# ---------------------
# Save standardized anomalies
# ---------------------
#standardized_output_path = '/share/data1/Students/jfields/finalproj/sst_anomalies_detrended_standardized_2020_2022.nc'
#standardized_anomalies.to_netcdf(standardized_output_path, format='NETCDF4_CLASSIC')

#print(f"Detrended and standardized anomalies saved to {standardized_output_path}")

output_path = '/share/data1/Students/jfields/finalproj/sst_anomalies_detrended_2020_2022.nc'
detrended_anomalies.to_netcdf(output_path, format='NETCDF4_CLASSIC')

print(f"Detrended anomalies saved to {output_path}")