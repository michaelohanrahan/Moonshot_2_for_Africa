import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point


def plot_annual_mean(pet, year, output_dir):
    annual_mean = pet.mean(dim='time')
    
    # Add cyclic point to prevent white line at 180 degrees longitude
    data_cyclic, lon_cyclic = add_cyclic_point(annual_mean.values, coord=annual_mean.x)
    
    # Create a new figure with a Robinson projection
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    
    # Set global extent and add coastlines
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    
    # Plot the data
    im = ax.pcolormesh(lon_cyclic, annual_mean.y, data_cyclic, 
                       transform=ccrs.PlateCarree(), 
                       cmap='viridis')
    
    # Add colorbar
    cbar = plt.colorbar(im, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('PET (mm/day)')
    
    # Set title
    plt.title(f'Annual Mean PET for {year}', fontsize=16)
    
    plt.show()
    # Save the figure
    # plot_path = Path(output_dir) / f'annual_mean_pet_{year}_robinson.png'
    # plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    # plt.close()

ds = xr.open_dataset(r"p:\moonshot2-casestudy\Wflow\africa\data\3-input\global_era5_pet\era5_daily_debruin_PET_daily_2024.nc")

# Count NaN values across all time steps
total_nan_count = np.isnan(ds.pet.values).sum()

print(f"Total number of NaN values across all time steps: {total_nan_count}")

# Calculate the percentage of NaN values across all time steps
total_elements = ds.pet.size
total_nan_percentage = (total_nan_count / total_elements) * 100

plot_annual_mean(ds.pet, 2024, "p:\\moonshot2-casestudy\\Wflow\\africa\\data\\3-input\\global_era5_pet\\")