import hydromt
import xarray as xr
import argparse
import traceback
from helper import setup_logging
from hydromt import DataCatalog
from icecream import ic
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
import yaml
import tempfile

def write_forcing(ds, fn_out, freq_out='Y', chunksize=1, decimals=2, time_units="days since 1900-01-01", **kwargs):
    """
    Write dataset to NetCDF files with a specified frequency (e.g., yearly, monthly, daily).
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
        if freq_out == 'Y':
            time_label = ds_gr["time"].dt.year[0].item()
        elif freq_out == 'M':
            time_label = f"{ds_gr['time'].dt.year[0].item()}-{ds_gr['time'].dt.month[0].item():02d}"
        elif freq_out == 'D':
            time_label = ds_gr["time"].dt.strftime("%Y-%m-%d").item()
        else:
            raise ValueError(f"Unsupported frequency: {freq_out}")

        fn_out_gr = fn_out.replace('*', str(time_label))
        if not os.path.isdir(os.path.dirname(fn_out_gr)):
            os.makedirs(os.path.dirname(fn_out_gr))
        
        delayed_obj = ds_gr.to_netcdf(fn_out_gr, encoding=encoding, mode="w", compute=False)
        print(f"{'#'*20}\nWRITING {freq_out.upper()}: {time_label} to: {fn_out_gr}\n{'#'*20}")
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

def plot_mean(args, pet=None, year=None, output_dir=None, month=None,):
    if month:
        #ic(pet)
        mean_data = pet.mean(dim='time')
        #ic(mean_data)
        title = f'Monthly Mean PET for {year}-{month:02d} using {args.method} method'
        filename = f'monthly_mean_pet_{year}_{month:02d}_{args.method}.png'
    else:
        mean_data = pet.mean(dim='time')
        title = f'Annual Mean PET for {year} using {args.method} method'
        filename = f'annual_mean_pet_{year}_{args.method}.png'
    
    # Add cyclic point to prevent white line at 180 degrees longitude
    data_cyclic, lon_cyclic = add_cyclic_point(mean_data.values, coord=mean_data.x)
    
    # Create a new figure with a Robinson projection
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    
    # Set global extent and add coastlines
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    
    # Plot the data
    im = ax.pcolormesh(lon_cyclic, mean_data.y, data_cyclic, 
                       transform=ccrs.PlateCarree(), 
                       cmap='viridis')
    
    # Add colorbar
    cbar = plt.colorbar(im, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('PET (mm/day)')
    
    # Set title
    plt.title(title, fontsize=16)
    
    # Save the figure
    plot_path = Path(output_dir) / filename
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
    temp_in = rename_dims(temp_in)
    # Convert pressure from Pa to hPa
    grid['press_msl'] = grid['press_msl'] / 100

    try:
        ic(temp_in)
        ic(grid['press_msl'])
        ic(grid['kin'])
        ic(grid['kout'])
        pet = pet_debruin(
            temp=temp_in,
            press=grid['press_msl'],
            k_in=grid['kin'],
            k_ext=grid['kout']
        )
        pet = pet.chunk(time=1)
        ic(pet)
        
        if np.isnan(pet.values).all():
            l.error(f"All PET values are NaN for the year. Debugging pet function...")
            raise ValueError("All PET values are NaN")
        
    except Exception as e:
        l.error(f"Error in PET calculation: {str(e)}")
        l.exception("Exception details:")
        pet = xr.full_like(grid, np.nan)

    return pet

def check_incomplete(grid, year, existing_file, month=None):
    if month:
        t0 = pd.Timestamp(f"{year}-{month:02d}-01")
        t1 = t0 + pd.offsets.MonthEnd()
        existing_file = xr.open_dataset(existing_file).sel(time=slice(t0, t1))
        #ic(existing_file)
        grid = grid.sel(time=slice(t0, t1))
        #ic(grid)
        if len(existing_file.time) < len(grid.time):
            return True
        else:
            return False
    else:
        if len(existing_file.time) < len(grid.time):
            return True
        else:
            return False

def process_single_month(grid, args, l, dc, outdir, year, month):
    # Convert month to integer if it's a string
    month = int(month)
    
    period_start = pd.Timestamp(f"{year}-{month:02d}-01")
    period_end = period_start + pd.offsets.MonthEnd()
    
    g_period = grid.sel(time=slice(period_start, period_end))
    
    if len(g_period.time) == 0:
        l.info(f"No data available for period {period_start} to {period_end}, skipping.")
        return
    
    # Correct the output file naming
    freq_out = 'hourly' if args.freq_data == 'H' else args.freq_data
    freq_out = 'daily' if freq_out == 'D' else freq_out
    
    output_file = Path(outdir) / f'{args.tpf}_{freq_out}_{args.method}_PET_{freq_out}_{year}_{month:02d}_{year}-{month:02d}.nc'
    
    if output_file.exists() and check_incomplete(g_period, year, output_file, month):
        l.info(f"Incomplete PET data for period {period_start} to {period_end}, overwrite set to True.")
        args.overwrite = True
    
    if output_file.exists() and not args.overwrite:
        l.warning(f"PET file for period {period_start} to {period_end} already exists, and is complete. Skipping. Use --overwrite to recalculate.")
        return
    
    l.info(f"Processing period {period_start} to {period_end}")
    pet = calc_pet(g_period, args, l, dc)
    l.info(f"PET calculated for period, time length: {len(pet.time.values)}")
    
    if isinstance(pet, xr.DataArray):
        # Plot monthly mean
        plot_mean(pet, year, Path(outdir) / args.method, month)
        l.info(f"Monthly mean plot created for {year}-{month:02d}")
        
        ds_out = create_output_dataset(pet, grid, args)
        
        l.info(f"Writing PET data for period {period_start} to {period_end}")
        write_forcing(ds_out, str(output_file), freq_out=args.freq_save, chunksize=1, decimals=2, time_units="days since 1900-01-01")
        l.info(f"Finished writing PET data for period {period_start} to {period_end}")
    else:
        l.error(f"Invalid pet data type: {type(pet)}")

def process_year(grid, args, l, dc, outdir, year):
    for month in range(1, 13):
        process_single_month(grid, args, l, dc, outdir, year, month)

def create_output_dataset(pet, grid, args, is_hourly=False):
    if is_hourly:
        timestep = 'Hourly'
        units = 'mm/hr'
    else:
        timestep = 'Daily'
        units = 'mm/day'
    return xr.Dataset({'pet': pet}, attrs={'description': f'Potential evapotranspiration calculated using {args.method} method with {args.tpf} temperature forcing',
                                           'meteo data source':f"{args.tpf}",
                                           'method': f"{args.method}",
                                           'dem': 'era5_orography',
                                           'timestep': timestep,
                                           'lapsecorrected': 'False',
                                           'presscorrected': 'False',
                                           'windaltitude': 'False',
                                           'units': units,
                                           'last updated': time.ctime(),
                                           **grid.attrs})

def adjust_yaml_root(yaml_path, drive):
    with open(yaml_path, 'r') as file:
        data = yaml.safe_load(file)
    
    # Adjust the root path in the meta section
    if 'meta' in data and 'root' in data['meta']:
        data['meta']['root'] = str(Path(drive) / 'wflow_global/hydromt')
    
    # Adjust paths in all top-level entries except 'meta'
    for key, value in data.items():
        if key != 'meta' and isinstance(value, dict) and 'path' in value:
            value['path'] = str(Path(data['meta']['root']) / value['path'].lstrip('/'))
    
    return data

def rename_dims(ds):
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'x'})
    if 'latitude' in ds.coords:
        ds=ds.rename({'latitude': 'y'})
    return ds

def main(args):
    os.chdir(Path(args.cwd))
    outdir = Path(args.cwd, 'data', '3-input', 'global_era5_pet').as_posix()
    os.makedirs(outdir, exist_ok=True)
    plot_dir = Path(outdir, args.method)
    os.makedirs(plot_dir, exist_ok=True)
    l = setup_logging('data/0-log', '00_build_PET.log')
    l.info("Building model assuming access to deltares_data catalog")
    
    drive = syscheck()
    yaml_path = Path(drive) / 'wflow_global/hydromt/deltares_data.yml'
    adjusted_yaml = adjust_yaml_root(yaml_path, drive)
    
    # Create a temporary file with the adjusted YAML
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as temp_file:
        yaml.dump(adjusted_yaml, temp_file)
        temp_yaml_path = temp_file.name
    
    try:
        dc = DataCatalog(temp_yaml_path)
        
        grid = dc.get_rasterdataset(args.tpf)
        grid = rename_dims(grid)  # Rename dimensions
        
        l.info(f"Calculating PET with method: {args.method}")
        
        # Determine the time step using a subset of the data
        time_subset = grid.time[:100]  # Take only the first 100 time steps
        time_diff = time_subset.diff('time').median()
        l.info(f"Median time difference: {time_diff}")
        
        is_hourly = time_diff <= np.timedelta64(1, 'h')
        
        if args.year is None:
            l.error("Year must be specified")
            return
        
        if args.month is None:
            l.info(f"Processing entire year {args.year}")
            process_year(grid, args, l, dc, outdir, args.year)
        else:
            l.info(f"Processing month {args.month} of year {args.year}")
            process_single_month(grid, args, l, dc, outdir, args.year, args.month)
        
        l.info("PET calculation completed for specified period.")
    finally:
        # Clean up the temporary file
        os.unlink(temp_yaml_path)

if __name__ == "__main__":
    try:
        drive = syscheck()
        parser = argparse.ArgumentParser(description='Build global PET with available data')
        parser.add_argument('--cwd', type=str, help='the current working directory', default=f'{drive}/moonshot2-casestudy/Wflow/africa')
        parser.add_argument('--tpf', type=str, help='the temperature forcing', default='era5_daily')
        parser.add_argument('--method', type=str, help='the method to use for PET calculation', default='debruin')
        parser.add_argument('--year', type=int, help='year to process', required=True)
        parser.add_argument('--month', type=str, help='month to process (optional, in format 02d)', default=None)
        parser.add_argument('--freq-data', type=str, help='data frequency', default='D')
        parser.add_argument('--freq-save', type=str, help='freq out', default='Y')
        parser.add_argument('--overwrite', action='store_true', help='overwrite existing PET files')
        args = parser.parse_args()
        main(args)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        traceback.print_exc()
    finally:
        print("Script finished.")