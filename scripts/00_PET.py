import hydromt
import xarray as xr
import argparse
import traceback
from helper import setup_logging
from hydromt import DataCatalog
import os
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar
import time
from helper import syscheck
from glob import glob
import numpy as np
import hydromt.workflows.forcing
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

def write_forcing(ds, fn_out, freq_out='Y', chunksize=1, decimals=2, time_units="days since 1900-01-01", **kwargs):
    """
    Write dataset to NetCDF files with a specified frequency (e.g., yearly).
    """
    if decimals is not None:
        ds = ds.round(decimals)

    chunksizes = (chunksize, ds.dims['y'], ds.dims['x'])
    encoding = {var: {"zlib": True, "dtype": "float32", "chunksizes": chunksizes} for var in ds.data_vars}

    for v in ["time", "x", "y"]:
        if v in ds.coords:
            ds[v].attrs.pop("_FillValue", None)
            encoding[v] = {"_FillValue": None}

    for label, ds_gr in ds.resample(time=freq_out):
        year = ds_gr["time"].dt.year[0].item()
        fn_out_gr = fn_out.replace('*', str(year))
        ##ic(fn_out_gr)
        if not os.path.isdir(os.path.dirname(fn_out_gr)):
            os.makedirs(os.path.dirname(fn_out_gr))
        
        delayed_obj = ds_gr.to_netcdf(fn_out_gr, encoding=encoding, mode="w", compute=False)
        print(f"{'#'*20}\nWRITING YEAR: {year} to: {fn_out_gr}\n{'#'*20}")
        with ProgressBar():
            delayed_obj.compute(**kwargs)

def vars(args):
        variables = ["temp"]
        if args.method == "debruin":
            variables += ["press_msl", "kin", "kout"]
        elif args.method == "makkink":
            variables += ["press_msl", "kin"]
        elif args.method == "penman-monteith_rh_simple":
            variables += ["temp_min", "temp_max", "rh", "kin"]
        elif args.method == "penman-monteith_tdew":
            variables += ["temp_min", "temp_max", "wind10_u", "wind10_v", "temp_dew", "kin", "press_msl"]
        return variables

def check_dates(outfile):
    files = glob(outfile)
    if not files:
        return None, None
    else:
        all_times = []
        for file in files:
            with xr.open_dataset(file) as ds:
                # Drop the spatial_ref coordinate if it exists
                if 'spatial_ref' in ds.coords:
                    ds = ds.drop('spatial_ref')
                if 'time' in ds.dims:
                    all_times.extend(ds.time.values)
        
        if not all_times:
            return None, None
        
        dates = pd.to_datetime(all_times)
        tmin = dates.min()
        tmax = dates.max()
        return tmin, tmax

def kelvin_to_celsius(temp_kelvin):
    return temp_kelvin - 273.15

def pet_debruin(temp, press, k_in, k_ext, timestep=86400, cp=1005.0, beta=20.0, Cs=110.0):
    esat = 6.112 * np.exp((17.67 * temp) / (temp + 243.5))
    slope = esat * (17.269 / (temp + 243.5)) * (1.0 - (temp / (temp + 243.5)))
    lam = (2.502 * 10**6) - (2250.0 * temp)
    gamma = (cp * press) / (0.622 * lam)
    ep_joule = (
        (slope / (slope + gamma))
        * (((1.0 - 0.23) * k_in) - (Cs * (k_in / (k_ext + 0.00001))))
    ) + beta
    ep_joule = xr.where(k_ext == 0.0, 0.0, ep_joule)
    pet = ((ep_joule / lam) * timestep).astype(np.float32)
    pet = xr.where(pet > 0.0, pet, 0.0)
    return pet

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
    
    # Save the figure
    plot_path = Path(output_dir) / f'annual_mean_pet_{year}_robinson.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

def calc_pet(grid, args, l, dc):
    DEM = dc.get_rasterdataset('era5_orography').isel(time=0, drop=True)
    l.info("Using ERA5 orography for DEM")
    assert DEM.raster.identical_grid(grid), "The grids are not identical"
    
    variables = vars(args)
    if 'u10' in grid.data_vars:
        grid = grid.rename({'u10': 'wind10_u', 'v10': 'wind10_v'})
    if 'd2m' in grid.data_vars:
        grid = grid.rename({'d2m': 'temp_dew'})

    temp_in = hydromt.workflows.forcing.temp(
        grid['temp'], dem_model=DEM, dem_forcing=DEM, lapse_correction=False, logger=l, freq=None)
    
    # Convert pressure from Pa to hPa
    grid['press_msl'] = grid['press_msl'] / 100

    try:
        pet = pet_debruin(
            temp=temp_in,
            press=grid['press_msl'],
            k_in=grid['kin'],
            k_ext=grid['kout']
        )
        
        if np.isnan(pet.values).all():
            l.error(f"All PET values are NaN for the year. Debugging pet function...")
        else:
            l.info(f"PET stats: min={pet.min().values.item():.2f}, max={pet.max().values.item():.2f}, mean={pet.mean().values.item():.2f}")
        
    except Exception as e:
        l.error(f"Error in PET calculation: {str(e)}")
        l.exception("Exception details:")
        pet = xr.full_like(grid, np.nan)

    return pet

def main(args):
    os.chdir(Path(args.cwd))
    outdir = Path(args.cwd, 'data', '3-input', 'global_era5_pet').as_posix()
    os.makedirs(outdir, exist_ok=True)
    plot_dir = Path(outdir, args.method)
    os.makedirs(plot_dir, exist_ok=True)
    outfile = Path(outdir, f'{args.tpf}_{args.method}_PET_daily_*.nc')
    
    l = setup_logging('data/0-log', '00_build_PET.log')
    l.info("Building model assuming access to deltares_data catalog")
    l.info(f"Writing output to {outfile}")
    
    drive = syscheck()
    dc = DataCatalog(f'{drive}/wflow_global/hydromt/deltares_data.yml')
    
    grid = dc.get_rasterdataset(args.tpf)
    
    l.info(f"Calculating PET with method: {args.method}")
    
    start_year = grid.time.dt.year.max().item()
    end_year = grid.time.dt.year.min().item()
    
    for year in range(start_year, end_year - 1, -1):
        year_start = pd.Timestamp(f"{year}-01-01")
        year_end = pd.Timestamp(f"{year}-12-31")
        
        g_year = grid.sel(time=slice(year_start, year_end))
        
        if len(g_year.time) == 0:
            l.info(f"No data available for year {year}, skipping.")
            continue
        
        output_file = Path(outdir) / f'{args.tpf}_{args.method}_PET_daily_{year}.nc'
        
        if output_file.exists() and not args.overwrite:
            l.warning(f"PET file for year {year} already exists. Skipping. Use --overwrite to recalculate.")
            continue
        
        l.info(f"Processing year {year}")
        pet = calc_pet(g_year, args, l, dc)
        l.info(f"PET calculated for year {year}, time length: {len(pet.time.values)}")
        
        ds_out = xr.Dataset({'pet': pet}, attrs={'description': f'Potential evapotranspiration calculated using {args.method} method with {args.tpf} temperature forcing',
                                                'meteo data source':f"{args.tpf}",
                                                'method': f"{args.method}",
                                                'dem': 'era5_orography',
                                                'timestep': 'daily',
                                                'lapsecorrected': 'False',
                                                'presscorrected': 'False',
                                                'windaltitude': 'False',
                                                'units': 'mm/day',
                                                'last updated': time.ctime(),
                                                **grid.attrs})
        
        ds_out = ds_out.rename({'latitude': 'y', 'longitude': 'x'})
        l.info(f"Writing PET data for year {year}")
        write_forcing(ds_out, str(output_file), freq_out='Y', chunksize=1, decimals=2, time_units="days since 1900-01-01T00:00:00")
        l.info(f"Finished writing PET data for year {year}")
        
        # Plot annual mean
        plot_annual_mean(pet, year, plot_dir)
        l.info(f"Annual mean plot created for year {year}")
    
    l.info("PET calculation completed for all available years.")

if __name__ == "__main__":
    try:
        drive = syscheck()
        parser = argparse.ArgumentParser(description='Build global PET with available data')
        parser.add_argument('--cwd', type=str, help='the current working directory', default=f'{drive}/moonshot2-casestudy/Wflow/africa')
        parser.add_argument('--tpf', type=str, help='the temperature forcing', default='era5_daily')
        parser.add_argument('--method', type=str, help='the method to use for PET calculation', default='debruin')
        parser.add_argument('--overwrite', action='store_true', help='overwrite existing PET files')
        args = parser.parse_args()
        main(args)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        traceback.print_exc()
    finally:
        print("Script finished.")