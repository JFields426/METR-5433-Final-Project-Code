import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

##################################################################################################################################
#REQUIRED USER INPUTS: 

#Line 15: Define file path for regional block ID dataset.
#Line 75: Define file path for output frequency map.
##################################################################################################################################

# Load the dataset
ds = xr.open_dataset('/file_path/regionblock.yyyymm1_yyyymm2.nc')

# Convert longitude to 0–360
ds = ds.assign_coords(longitude=(ds.longitude % 360))

# Subset to 35N–80N, 290E–360E
subset = ds.sel(latitude=slice(35, 80), longitude=slice(290, 360))

# Extract object_id
object_id = subset['object_id']

# Total number of time steps (days)
total_days = object_id.sizes['time']

# Count non-zero occurrences over time
nonzero_counts = (object_id != 0).sum(dim='time')

# Calculate percentage
percentage_nonzero = (nonzero_counts / total_days) * 100

# Convert longitude back to -180–180 for plotting labels
plot_longitudes = ((percentage_nonzero.longitude + 180) % 360) - 180
plot_lats = percentage_nonzero.latitude

# Create meshgrid
lon2d, lat2d = np.meshgrid(plot_longitudes, plot_lats)

# Plot with Cartopy
plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot data with white-to-blue colormap
cmap = plt.cm.plasma
mesh = ax.pcolormesh(
    lon2d,
    lat2d,
    percentage_nonzero,
    cmap=cmap,
    shading='auto',
    transform=ccrs.PlateCarree()
)

# Add coastlines and land features
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='lightgray')

# Set extent to the region of interest (-70° to 0°, 35N–80N)
ax.set_extent([-70, 0, 35, 80], crs=ccrs.PlateCarree())

# Add colorbar
cbar = plt.colorbar(mesh, orientation='vertical', pad=0.05, shrink=0.7)
cbar.set_label('Frequency of Blocked Grid Points (%)')

# Add title and labels
plt.title('Blocking Frequency in the Northern Atlantic Region\n(35°N-80°N, 70°W-0°E) (1950–2022)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Save plot
plt.savefig('frequency_map.yyyymm1_yyyymm2.png', dpi=300, bbox_inches='tight')
plt.close()
