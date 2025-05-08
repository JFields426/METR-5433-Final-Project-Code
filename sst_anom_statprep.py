import xarray as xr

#########################################################################
#REQUIRED USER INPUTS: 

#Line 10: Define file path and name for the input file
#Line 38: Define file path and name for the output file
#########################################################################

anom_input_path = '/file_path/sst_anomalies.yyyymm1_yyyymm2.nc'
anomalies = xr.open_dataset(anom_input_path)['sst']

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

# Standardize Anomalies
mean = detrended_anomalies.mean(dim='time', skipna=True)
std = detrended_anomalies.std(dim='time', skipna=True)

standardized_anomalies = (detrended_anomalies - mean) / std

# Name the variable explicitly
standardized_anomalies.name = 'sst_anom'

# Save standardized anomalies
standardized_output_path = '/file_path/sst_anomalies_detrended_standardized.yyymm1_yyyymm2.nc'
standardized_anomalies.to_netcdf(standardized_output_path, format='NETCDF4_CLASSIC')

print(f"Detrended and standardized anomalies saved to {standardized_output_path}")
