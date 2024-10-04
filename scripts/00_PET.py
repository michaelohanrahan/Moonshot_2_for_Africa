import hydromt
import xarray as xr
import argparse
import traceback
from helper import setup_logging
from hydromt_wflow import WflowModel as WM
from hydromt import DataCatalog
import os
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar
import time
from helper import syscheck
from glob import glob
from tqdm import tqdm
from icecream import ic

def write_forcing(ds, fn_out, freq_out='Y', chunksize=1, decimals=2, time_units="days since 1900-01-01T00:00:00", **kwargs):
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

        if not os.path.isdir(os.path.dirname(fn_out_gr)):
            os.makedirs(os.path.dirname(fn_out_gr))
        
        delayed_obj = ds_gr.to_netcdf(fn_out_gr, encoding=encoding, mode="w", compute=False)
        
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
    if glob(outfile)==[]:
        return None, None
    else:
        ds = xr.open_mfdataset(outfile)
        dates = pd.to_datetime(ds.time.values)
        tmin = dates.min()
        tmax = dates.max()
        return tmin, tmax

def calc_pet(grid, args, l, dc):
    DEM = dc.get_rasterdataset(
        'era5_orography'
    ).isel(time=0, drop=True)
    
    l.warning("Using ERA5 orography for DEM")

    assert DEM.raster.identical_grid(grid) == True, "The grids are not identical"
    
    variables = vars(args)
    if 'u10' in grid.data_vars:
        grid = grid.rename({'u10': 'wind10_u', 'v10': 'wind10_v'})
    if 'd2m' in grid.data_vars:
        grid = grid.rename({'d2m': 'temp_dew'})

    temp_in = hydromt.workflows.forcing.temp(
        grid['temp'], dem_model=DEM, dem_forcing=DEM, lapse_correction=False, logger=l, freq=None)
    if "penman-monteith" in args.method:
        l.info("Calculating max and min temperature")
        temp_max_in = hydromt.workflows.forcing.temp(
            grid['temp_max'], dem_model=DEM, dem_forcing=DEM, lapse_correction=False, logger=l, freq=None).rename('temp_max')
        temp_min_in = hydromt.workflows.forcing.temp(
            grid['temp_min'], dem_model=DEM, dem_forcing=DEM, lapse_correction=False, logger=l, freq=None).rename('temp_min')
        temp_in = xr.merge([temp_in, temp_max_in, temp_min_in])
    year = grid.time.dt.year[0].item()
    das = []
    for t in tqdm(grid.time.values, desc=f'Calculating PET {year}', unit='day', position=0, total=len(grid.time.values), colour='green'):
        day = hydromt.workflows.forcing.pet(
            grid[variables[1:]].sel(time=slice(t,t), drop=False), temp=temp_in.sel(time=slice(t,t), drop=False), dem_model=DEM, method=args.method, press_correction=False,
            wind_correction=False, wind_altitude=False, freq='D')
        das.append(day)
    pet = xr.concat(das, dim='time')
    l.info(f"Finished calculating PET for {year}")
    return pet

def main(args):
    os.chdir(Path(args.cwd))
    # Read configuration and set up output directory
    outdir = Path(args.cwd, 'data', '4-output', 'PET_global').as_posix()
    os.makedirs(outdir, exist_ok=True)
    outfile = Path(outdir, f'{args.tpf}_{args.method}_PET_daily_*.nc')
    
    # Check existing data dates
    tmin, tmax = check_dates(str(outfile))
    
    # Set up logging
    l = setup_logging('data/0-log', '00_build_PET.log')
    l.info("Building model assuming access to deltares_data catalog")
    l.info(f"Writing output to {outfile}")
    
    # Initialize data catalog
    drive = syscheck()
    dc = DataCatalog(f'{drive}/wflow_global/hydromt/deltares_data.yml')
    
    # Load input data based on existing dates or from the beginning
    if tmin is not None:
        l.info(f"Calculating PET from {tmin} to {tmax}")
        grid = dc.get_rasterdataset(
            args.tpf
        ).sel(time=slice(tmin,tmax)).chunk({"time": 1, "latitude": -1, "longitude": -1})
    else:
        l.info("Calculating PET from the beginning of the dataset")
        grid = dc.get_rasterdataset(
            args.tpf
        ).chunk({"time": 1, "latitude": -1, "longitude": -1})
        
    l.info(f"calculating PET with method: {args.method}")
    
    # Process data year by year
    year_range = pd.date_range(start=grid.time.values[0], end=grid.time.values[-1], freq='Y')
    
    for year in year_range:
        # Calculate PET for the current year
        g_year = grid.sel(time=slice(year, year+pd.Timedelta(days=365)))
        pet = calc_pet(g_year, args, l, dc)
        
        # Prepare output dataset with metadata
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
        
        # Rename coordinates and write output
        ds_out = ds_out.rename({'latitude': 'y', 'longitude': 'x'})
        write_forcing(ds_out, str(outfile), freq_out='Y', chunksize=1, decimals=2, time_units="days since 1900-01-01T00:00:00")
    outdir = Path(args.cwd, 'data', '4-output', 'PET_global').as_posix()
    os.makedirs(outdir, exist_ok=True)
    outfile = Path(outdir, f'{args.tpf}_{args.method}_PET_daily_*.nc')
    
    #TODO: Make it so that we append the new data to the existing file
    tmin, tmax = check_dates(str(outfile))
    
    l = setup_logging('data/0-log', '00_build_PET.log')
    l.info("Building model assuming access to deltares_data catalog")
    l.info(f"Writing output to {outfile}")
    
    drive = syscheck()
    
    dc = DataCatalog(f'{drive}/wflow_global/hydromt/deltares_data.yml')
    
    if tmin is not None:
        l.info(f"Calculating PET from {tmin} to {tmax}")
        grid = dc.get_rasterdataset(
            args.tpf
        ).sel(time=slice(tmin,tmax)).chunk({"time": 1, "latitude": -1, "longitude": -1})
        
    else:
        l.info("Calculating PET from the beginning of the dataset")
        grid = dc.get_rasterdataset(
            args.tpf
        ).chunk({"time": 1, "latitude": -1, "longitude": -1})
        
    l.info(f"calculating PET with method: {args.method}")
    l.info(f"Calculating PET for {grid.time.values[0]} to {grid.time.values[-1]}")
    l.info(f"Using {args.tpf} as the temperature forcing")
    l.info(f"Using {args.method} as the PET method")
    l.info(f"Using {args.dem} as the DEM")  
    l.info(f"Using {args.timestep} as the timestep")
    l.info(f"Lapse correction set to: {args.lapsecorrected}")
    l.info(f"Pressure correction set to: {args.presscorrected}")
    l.info(f"Wind altitude correction set to: {args.windaltitude}")
    l.info(f"Using {args.lastupdated} as the lastupdated")
    
    year_range = pd.date_range(start=grid.time.values[0], end=grid.time.values[-1], freq='Y')
    
    for year in year_range:
        g_year = grid.sel(time=slice(year, year+pd.Timedelta(days=365)))
        pet = calc_pet(g_year, args, l, dc)
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
        write_forcing(ds_out, str(outfile), freq_out='Y', chunksize=1, decimals=2, time_units="days since 1900-01-01T00:00:00")
    
    #write dataflags in same folder as readme
    #if files are not full years, write a dataflag file
    files = glob(str(outfile))
    
    for file in files:
        ds = xr.open_dataset(file)
        if ds.time.size < 365:
            with open(os.path.join(outdir, 'dataflags.txt'), 'a') as f:
                f.write(f'{file} does not contain a full 365 days\n')
        ds.close()

if __name__ == "__main__":
    try:
        drive = syscheck()
        parser = argparse.ArgumentParser(description='Build global PET with available data')
        parser.add_argument('--cwd', type=str, help='the current working directory', default=f'{drive}/moonshot2-casestudy/Wflow/africa')
        parser.add_argument('--tpf', type=str, help='the temperature forcing', default='era5_daily_zarr')
        parser.add_argument('--method', type=str, help='the method to use for PET calculation', default='penman-monteith_tdew')
        args = parser.parse_args()

        main(args)

    except SystemExit:
        traceback.print_exc()
