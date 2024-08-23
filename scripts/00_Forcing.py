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

def main(args):
    
    
    os.chdir(Path(r'p:\moonshot2-casestudy\Wflow\africa'))
    l = setup_logging('data/0-log', '00_forcing.log')
    l.info("Building model assuming access to deltares_data catalog")

    root = f'src/3-model/wflow_build/{args.clusterid}'
    config = f"wflow_sbm.toml"

    mod = WM(root=root,
             mode='r',
             config_fn=config,
             data_libs=['p:/wflow_global/hydromt/deltares_data.yml'],
             logger=l)

    bbox = mod.region.geometry
    bbox_b = bbox.buffer(0.5).bounds.values[0]

    dc = DataCatalog('p:/wflow_global/hydromt/deltares_data.yml')

    xmin, ymin, xmax, ymax = bbox_b
    
    tmax = args.tmax
    

    if ymin > ymax:
        # Case where ymin is greater than ymax
        grid = dc.get_rasterdataset(
            args.tpf
        ).sel(
            longitude=slice(xmin, xmax), latitude=slice(ymin, ymax)
        ).chunk({"time": 1})
        if tmax is None:
            tmax = pd.to_datetime(grid.time.values[-1])
        grid = grid.sel(time=slice(args.tmin, tmax))
        DEM = dc.get_rasterdataset(
            'era5_orography'
        ).sel(
            longitude=slice(xmin, xmax), latitude=slice(ymin, ymax)
        ).isel(time=0, drop=True)
    else:
        # Case where ymax is greater than ymin
        grid = dc.get_rasterdataset(
            args.tpf
        ).sel(
            longitude=slice(xmin, xmax), latitude=slice(ymax, ymin)
        ).chunk({"time": 1})
        if tmax is None:
            tmax = pd.to_datetime(grid.time.values[-1])
        grid = grid.sel(time=slice(args.tmin, tmax))
        DEM = dc.get_rasterdataset(
            'era5_orography'
        ).sel(
            longitude=slice(xmin, xmax), latitude=slice(ymax, ymin)
        ).isel(time=0, drop=True)


    assert DEM.raster.identical_grid(grid) == True, "The grids are not identical"

    variables = ["temp"]
    if args.method == "debruin":
        variables += ["press_msl", "kin", "kout"]
    elif args.method == "makkink":
        variables += ["press_msl", "kin"]
    elif args.method == "penman-monteith_rh_simple":
        variables += ["temp_min", "temp_max", "rh", "kin"]
    elif args.method == "penman-monteith_tdew":
        variables += ["temp_min", "temp_max", "wind10_u", "wind10_v", "temp_dew", "kin", "press_msl"]

    if 'u10' in grid.data_vars:
        grid = grid.rename({'u10': 'wind10_u', 'v10': 'wind10_v'})
    if 'd2m' in grid.data_vars:
        grid = grid.rename({'d2m': 'temp_dew'})

    P = grid['precip']
    precip = hydromt.workflows.forcing.precip(P, da_like=P, freq='D', logger=l)

    temp_in = hydromt.workflows.forcing.temp(
        grid['temp'], dem_model=DEM, dem_forcing=DEM, lapse_correction=False, logger=l, freq=None)

    if "penman-monteith" in args.method:
        l.info("Calculating max and min temperature")
        temp_max_in = hydromt.workflows.forcing.temp(
            grid['temp_max'], dem_model=DEM, dem_forcing=DEM, lapse_correction=False, logger=l, freq=None).rename('temp_max')
        temp_min_in = hydromt.workflows.forcing.temp(
            grid['temp_min'], dem_model=DEM, dem_forcing=DEM, lapse_correction=False, logger=l, freq=None).rename('temp_min')
        temp_in = xr.merge([temp_in, temp_max_in, temp_min_in])

    pet = hydromt.workflows.forcing.pet(
        grid[variables[1:]], temp=temp_in, dem_model=DEM, method=args.method, press_correction=False,
        wind_correction=False, wind_altitude=False, freq='D', logger=l)

    if "penman-monteith" in args.method:
        temp_in = temp_in.drop(['temp_max', 'temp_min'])

    temp_out = hydromt.workflows.forcing.resample_time(
        temp_in, freq='D', upsampling="bfill", downsampling="mean", label="right", closed="right", conserve_mass=False, logger=l)

    ds_out = xr.merge([temp_out, pet, precip]).chunk({"time": 1})
    ds_out = ds_out.rename({'latitude': 'y', 'longitude': 'x'})
    os.makedirs(os.path.join(mod.root, "inmaps"), exist_ok=True)
    outdir = f'inmaps/inmaps_pet_prec_{args.tpf}_{args.method}_daily_*.nc'

    write_forcing(ds_out, os.path.join(mod.root, outdir), freq_out='Y', chunksize=1, decimals=2, time_units="days since 1900-01-01T00:00:00")

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description='Build Forcing for a Bounding Box')
        parser.add_argument('--tpf', type=str, help='the temperature forcing', default='era5_daily_zarr')
        parser.add_argument('--tpr', type=str, help='the precipitation forcing', default='era5_daily_zarr')
        parser.add_argument('--method', type=str, help='the method to use for PET calculation', default='penman-monteith_tdew')
        parser.add_argument('--clusterid', type=str, help='the cluster id', default='default_model')
        parser.add_argument('--tmin', type=str, help='the start time for the PET forcing %Y-%m-%d', default='1980')
        parser.add_argument('--tmax', type=str, help='the end time for the PET forcing %Y-%m-%d', default=None)
        args = parser.parse_args()

        main(args)

    except SystemExit:
        # This block is for running in environments like Jupyter
        traceback.print_exc()
