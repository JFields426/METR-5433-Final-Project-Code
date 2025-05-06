import xarray as xr

#########################################################################
#REQUIRED USER INPUTS: 

#Line 10: Define file path and name for the input file
#Line 31: Define file path and name for the output file
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

detrended_anomalies.name = 'sst_anom'

#Save dataset
output_path = '/file_path/sst_anomalies_detrended.yyyymm1_yyyymm2.nc'
detrended_anomalies.to_netcdf(output_path, format='NETCDF4_CLASSIC')

print(f"Detrended anomalies saved to {output_path}")
