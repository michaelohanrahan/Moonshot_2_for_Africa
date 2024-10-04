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
from helper import syscheck as sc
import matplotlib.pyplot as plt
import pyet
import numpy as np

#FUNCTIONS
def get_variables(args):
    variables = ["temp"]
    if args.method == "debruin":
        variables += ["press_msl", 
                    "kin",
                    "kout"]
    elif args.method == "makkink":
        variables += ["press_msl", 
                    "kin"]
    elif args.method == "penman-monteith_rh_simple":
        variables += ["temp_min", 
                    "temp_max", 
                    "rh", 
                    "kin"]
    elif args.method == "penman-monteith_tdew":
            variables += ["temp_max", 
                    "temp_min", 
                    "wind10_u",
                    "wind10_v", 
                    "temp_dew", 
                    "kin", 
                    "press_msl"]
    else:
        raise ValueError(f"Method {args.method} not recognized, must be one of: debruin, makkink, penman-monteith_rh_simple, penman-monteith_tdew")
    return variables
def write_forcing(l,
                  ds:xr.Dataset, 
                  fn_out:str, 
                  freq_out:str='Y', 
                  chunksize:int=1, 
                  decimals:int=2, 
                  time_units:str="days since 1900-01-01T00:00:00", 
                  **kwargs):
    """
    Write dataset to NetCDF files with a specified frequency (e.g., yearly).
    """
    l.info(f"Writing PET to {fn_out}")
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
        # l.info(f"Writing forcing for year {year}")
        # l.info(f"ds_gr.time.min().values: {ds_gr.time.min().values}")
        # l.info(f"ds_gr.time.max().values: {ds_gr.time.max().values}")
        fn_out_gr = fn_out.replace('*', str(year))

        if not os.path.isdir(os.path.dirname(fn_out_gr)):
            os.makedirs(os.path.dirname(fn_out_gr))
        
        delayed_obj = ds_gr.to_netcdf(fn_out_gr, encoding=encoding, mode="w", compute=False)
        
        with ProgressBar():
            delayed_obj.compute(**kwargs)
def model_bounds(l, args, datalib):
    ymin, ymax, xmin, xmax = None, None, None, None
    mod = None
    if args.clusterid is not None:
        l.info(f"Getting model bounds for cluster {args.clusterid}")
        root = f'src/3-model/wflow_build/{args.clusterid}'
        config = f"wflow_sbm.toml"

        mod = WM(root=root,
                mode='r',
                config_fn=config,
                data_libs=[datalib],
                logger=l)

        bbox = mod.region.geometry
        bbox_b = bbox.buffer(0.5).bounds.values[0]
        
        xmin, ymin, xmax, ymax = bbox_b
    elif args.bbox is not None:
        xmin, ymin, xmax, ymax = args.bbox.strip('[]').split(',')
        xmin, ymin, xmax, ymax = float(xmin)-1, float(ymin)-1, float(xmax)+1, float(ymax)+1
        
    # return 0 , 5 , 0 , 5 , mod
    return ymin, ymax, xmin, xmax, mod 
def slice_grid(grid, DEM, tmin, tmax, xmin, xmax, ymin, ymax):
    
    if tmax is None:
        tmax = pd.to_datetime(grid.time.values[-1])
        l.info(f"No tmax supplied, using full availability up to {tmax}")

    grid = grid.sel(time=slice(tmin, tmax))

    if ymin is not None and ymax is not None:
        ymin = grid.sel(latitude=ymin, method='nearest').latitude.values
        ymax = grid.sel(latitude=ymax, method='nearest').latitude.values
        xmin = grid.sel(longitude=xmin, method='nearest').longitude.values
        xmax = grid.sel(longitude=xmax, method='nearest').longitude.values
        ic(xmin, xmax, ymin, ymax)
        ic(grid.longitude.max().values, grid.longitude.min().values)
        ic(grid.latitude.max().values, grid.latitude.min().values)
        grid = grid.sel(
            longitude=slice(xmin, xmax), latitude=slice(ymax, ymin)
        )
        DEM = DEM.sel(
            longitude=slice(xmin, xmax), latitude=slice(ymax, ymin)
        )
    else:
        l.info(f"Computing global PET for {(tmin, tmax)}")
    
    return grid, DEM
# def slice_grid(grid, DEM, tmin, tmax, xmin, xmax, ymin, ymax):
#     if tmax is None:
#         tmax = pd.to_datetime(grid.time.values[-1])
#         l.info(f"No tmax supplied, using full availability up to {tmax}")

#     grid = grid.sel(time=slice(tmin, tmax))

#     if ymin is not None and ymax is not None:
#         # Convert longitudes to ±180 system
#         grid = grid.assign_coords(longitude=(((grid.longitude + 180) % 360) - 180))
#         DEM = DEM.assign_coords(longitude=(((DEM.longitude + 180) % 360) - 180))
#         # Ensure longitudes are sorted
#         grid = grid.sortby('longitude')
#         DEM = DEM.sortby('longitude')
        
#         ymin = grid.sel(latitude=ymin, method='nearest').latitude.values
#         ymax = grid.sel(latitude=ymax, method='nearest').latitude.values
#         xmin = grid.sel(longitude=xmin, method='nearest').longitude.values
#         xmax = grid.sel(longitude=xmax, method='nearest').longitude.values
        
#        # Select the region
#         if xmin <= xmax:
#             grid = grid.sel(longitude=slice(xmin, xmax), latitude=slice(ymax, ymin))
#             DEM = DEM.sel(longitude=slice(xmin, xmax), latitude=slice(ymax, ymin))
#         else:
#             # Handle case where selection crosses the antimeridian
#             grid = xr.concat([grid.sel(longitude=slice(xmin, 180)), 
#                               grid.sel(longitude=slice(-180, xmax))], dim='longitude')
#             DEM = xr.concat([DEM.sel(longitude=slice(xmin, 180)), 
#                              DEM.sel(longitude=slice(-180, xmax))], dim='longitude')
        
#         grid = grid.sel(latitude=slice(ymax, ymin))
#         DEM = DEM.sel(latitude=slice(ymax, ymin))

#         # Apply nearest neighbor selection if needed
#         grid = grid.sel(longitude=grid.longitude, latitude=grid.latitude, method='nearest')
#         DEM = DEM.sel(longitude=DEM.longitude, latitude=DEM.latitude, method='nearest')

#         # Ensure the grid is regular
#         lon_step = np.diff(grid.longitude).mean()
#         lat_step = np.diff(grid.latitude).mean()
#         new_lon = np.arange(grid.longitude.min().item(), grid.longitude.max().item() + lon_step, lon_step)
#         new_lat = np.arange(grid.latitude.min().item(), grid.latitude.max().item() + lat_step, lat_step)
#         grid = grid.interp(longitude=new_lon, latitude=new_lat)
#         DEM = DEM.interp(longitude=new_lon, latitude=new_lat)

#     else:
#         l.info(f"Computing global PET for {(tmin, tmax)}")
    
#     # Add raster attributes
#     grid.attrs['transform'] = (lon_step, 0.0, grid.longitude.min().item(),
#                                0.0, -lat_step, grid.latitude.max().item())
#     grid.attrs['res'] = (lon_step, lat_step)
#     grid.attrs['crs'] = 'EPSG:4326'

#     DEM.attrs['transform'] = grid.attrs['transform']
#     DEM.attrs['res'] = grid.attrs['res']
#     DEM.attrs['crs'] = grid.attrs['crs']

#     fig, ax = plt.subplots()
#     ax.set_title(f"debugDEM")
#     DEM.plot(ax=ax)
#     plt.savefig('debugDEM.png')
#     plt.close()

#     return grid, DEM

def prep_vars(l, grid, DEM, args):
    if 'u10' in grid.data_vars:
        grid = grid.rename({'u10': 'wind10_u', 
                            'v10': 'wind10_v'})
        
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
    
    temp_in = temp_ktoc(temp_in)
    grid = temp_ktoc(grid)
    
    if 'rh' in args.method.split('_'):
        l.info("Calculating relative humidity")
        grid['rh'] = 5 * grid['temp_dew'] + 100 - 5 * grid['temp']
    #check units and convert to celsius if necessary
    
    return temp_in, grid
def calc_pet(
    ds: xr.Dataset,
    temp: xr.DataArray,
    dem_model: xr.DataArray,
    method: str = "debruin",
    press_correction: bool = False,
    wind_correction: bool = True,
    wind_altitude: float = 10,
    reproj_method: str = "nearest_index",
    freq: str = None,
    resample_kwargs: dict = None,
    logger=None,
) -> xr.DataArray:
    """Determine reference evapotranspiration.

    (lazy reprojection on model grid and resampling of time dimension to frequency).

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with climate variables: pressure [hPa], global radiation [W m-2],
        TOA incident solar radiation [W m-2], wind [m s-1]

        * Required variables: {"temp", "press" or "press_msl", "kin"}
        * additional variables for debruin: {"kout"}
        * additional variables for penman-monteith_rh_simple:
            {"temp_min", "temp_max", "wind" or "wind_u"+"wind_v", "rh"}
        * additional variables for penman-monteith_tdew:
            {"temp_min", "temp_max", "wind" or "wind_u"+"wind_v", "temp_dew"}
    temp : xarray.DataArray
        DataArray with temperature on model grid resolution [°C]
    dem_model : xarray.DataArray
        DataArray of the target resolution and projection, contains elevation
        data
    method : {'debruin', 'makkink', "penman-monteith_rh_simple", "penman-monteith_tdew"}
        Potential evapotranspiration method.
        if penman-monteith is used, requires the installation of the pyet package.
    reproj_method: str, optional
        Method for spatital reprojection of precip, by default 'nearest_index'
    press_correction : bool, default False
        If True pressure is corrected, based on elevation data of `dem_model`
    wind_altitude: float, optional
        ALtitude of wind speed data. By default 10m.
    wind_correction: str, optional
        If True wind speed is re-calculated to wind speed at 2 meters using
        original `wind_altitude`.
    freq: str, Timedelta, default None
        output frequency of timedimension
    resample_kwargs:
        Additional key-word arguments (e.g. label, closed) for time resampling method
    logger:
        The logger to use.

    Returns
    -------
    pet_out : xarray.DataArray (lazy)
        reference evapotranspiration
    """
    resample_kwargs=resample_kwargs or {}
    # resample input to model grid
    ds = ds.raster.reproject_like(dem_model, method=reproj_method)
    # logger.info(f"debug: match all lat ds.latitude.values: {np.all(ds.latitude.values == dem_model.latitude.values)}")
    # logger.info(f"debug: match all lon ds.longitude.values: {np.all(ds.longitude.values == dem_model.longitude.values)}")
    # Process bands like 'pressure' and 'wind'
    if press_correction:
        ds["press"] = hydromt.workflows.forcing.press(
            ds["press_msl"],
            dem_model,
            lapse_correction=press_correction,
            freq=None,  # do not change freq of press, put pet_out later
            reproj_method=reproj_method,
        )
    else:
        if "press_msl" in ds:
            ds = ds.rename({"press_msl": "press"})
            # logger.info(f"debug: min press: {ds['press'].min().values}, max press: {ds['press'].max().values}, mean press: {ds['press'].mean().values}")
        elif pyet.__version__ is not None:
            # calculate pressure from elevation [kPa]
            ds["press"] = pyet.calc_press(dem_model)
            # convert to hPa to be consistent with press function calculation:
            ds["press"] = ds["press"] * 10
        else:
            raise ModuleNotFoundError(
                "If 'press' is supplied and 'press_correction' is not used,"
                + " the pyet package must be installed."
            )

    timestep = pd.to_timedelta(ds.time.values[1] - ds.time.values[0]).total_seconds()
    
    if method == "debruin":
        pet_out = hydromt.workflows.forcing.pet_debruin(
            temp,
            ds["press"],
            ds["kin"],
            ds["kout"],
            timestep=timestep,
        )
        # logger.info(f"debug debruin: min pet: {pet_out.min().values}, max pet: {pet_out.max().values}, mean pet: {pet_out.mean().values}")
        
    elif method == "makkink":
        pet_out = hydromt.workflows.forcing.pet_makkink(temp, ds["press"], ds["kin"], timestep=timestep)
        # logger.info(f"debug makkink: min pet: {pet_out.min().values}, max pet: {pet_out.max().values}, mean pet: {pet_out.mean().values}")
        
    elif "penman-monteith" in method:
        logger.info("Calculating Penman-Monteith ref evaporation")
        # Add wind
        # compute wind from u and v components at 10m (for era5)
        if ("wind10_u" in ds.data_vars) & ("wind10_v" in ds.data_vars):
            ds["wind"] = hydromt.workflows.forcing.wind(
                da_model=dem_model,
                wind_u=ds["wind10_u"],
                wind_v=ds["wind10_v"],
                altitude=wind_altitude,
                altitude_correction=wind_correction,
            )
            # logger.info(f"debug wind: min wind: {ds['wind'].min().values}, max wind: {ds['wind'].max().values}, mean wind: {ds['wind'].mean().values}")
        else:
            ds["wind"] = hydromt.workflows.forcing.wind(
                da_model=dem_model,
                wind=ds["wind"],
                altitude=wind_altitude,
                altitude_correction=wind_correction,
            )
            # logger.info(f"debug wind: min wind: {ds['wind'].min().values}, max wind: {ds['wind'].max().values}, mean wind: {ds['wind'].mean().values}")
        if method == "penman-monteith_rh_simple":
            pet_out = hydromt.workflows.forcing.pm_fao56(
                temp["temp"],
                temp["temp_max"],
                temp["temp_min"],
                ds["press"],
                ds["kin"],
                ds["wind"],
                ds["rh"],
                dem_model,
                "rh",
            )
            # logger.info(f"debug pm_fao56: min pet: {pet_out.min().values}, max pet: {pet_out.max().values}, mean pet: {pet_out.mean().values}")
        elif method == "penman-monteith_tdew":
            pet_out = hydromt.workflows.forcing.pm_fao56(
                temp["temp"],
                temp["temp_max"],
                temp["temp_min"],
                ds["press"],
                ds["kin"],
                ds["wind"],
                ds["temp_dew"],
                dem_model,
                "temp_dew",
            )
            # logger.info(f"debug pm_fao56: min pet: {pet_out.min().values}, max pet: {pet_out.max().values}, mean pet: {pet_out.mean().values}")
        else:
            methods = [
                "debruin",
                "makking",
                "penman-monteith_rh_simple",
                "penman-monteith_tdew",
            ]
            raise ValueError(f"Unknown pet method, select from {methods}")

    # resample in time
    pet_out.name = "pet"
    pet_out.attrs.update(unit="mm")
    logger.info(f"debug pet_out: min pet: {pet_out.min().values}, max pet: {pet_out.max().values}, mean pet: {pet_out.mean().values}")
    if freq is not None:
        resample_kwargs.update(upsampling="bfill", downsampling="sum", logger=logger)
        pet_out = hydromt.workflows.forcing.resample_time(pet_out, freq, conserve_mass=True, **resample_kwargs)
    
    return pet_out

def temp_ktoc(temp):
    if isinstance(temp, xr.Dataset):
        for var in temp.data_vars:
            if 'temp' in var.split('_'):
                l.info(f"{var} is in K and all methods require C.. conversion applied")
                data_point = temp[var].isel(time=0, latitude=0, longitude=0).values
                if data_point >= 273.15:
                    temp[var] = temp[var] - 273.15
                    temp[var].attrs['units'] = 'C'
    
    elif isinstance(temp, xr.DataArray):
        if temp.sel(time=0, latitude=0, longitude=0).values >= 273.15:
            temp = temp - 273.15
            temp.attrs['units'] = 'C'
    return temp

#MAIN
def main(l, args):
    DRIVE=sc()
    os.chdir(Path(f'{DRIVE}/moonshot2-casestudy/Wflow/africa'))
    l.info("Building model assuming access to deltares_data catalog")
    if ':' in DRIVE:
        datalib = 'p:/wflow_global/hydromt/deltares_data.yml'
        dc = DataCatalog(datalib)
    else:
        datalib = Path(os.getcwd(), 'data/1-external/deltares_data_linux.yml').as_posix()
        dc = DataCatalog(datalib)
        
    ymin, ymax, xmin, xmax, mod = model_bounds(l, args, datalib)

    tmin = pd.Timestamp(args.tmin)
    tmax = pd.Timestamp(args.tmax)
    if tmin == tmax:
        stamp = f'{tmax.year}-12-31 00:00:00'
        l.warning(f" *** tmin and tmax are in the same year, setting tmax to {stamp} *** \n *** to avoid pass a string 'yyyy-mm-dd' ***")
        tmax = pd.Timestamp(stamp)
    
    l.info(f"Extracting PET for {tmin} to {tmax}")
    
    if tmax is None:
        l.info(f"No tmax provided, using most recent available year")
    
    grid = dc.get_rasterdataset(args.tpf).chunk({"time": 1})	#, "latitude": 100, "longitude": 100})

    DEM = dc.get_rasterdataset(
        'era5_orography'
    ).isel(time=0, drop=True)#.chunk({"latitude": 100, "longitude": 100})

    if ymin is None:
        ymin = DEM.latitude.values.min()
    if ymax is None:
        ymax = DEM.latitude.values.max()
    if xmin is None:
        xmin = DEM.longitude.values.min()
    if xmax is None:
        xmax = DEM.longitude.values.max()
    
    grid, DEM = slice_grid(grid, DEM, tmin, tmax, xmin, xmax, ymin, ymax)

    assert DEM.raster.identical_grid(grid) == True, "The forcing and DEM grids are not identical"
    
    variables = get_variables(args)
    
    temp_in, grid = prep_vars(l, grid, DEM, args)
    
    # for var in temp_in.data_vars:
    #     l.info(f"debug {var}: min {temp_in[var].min().values}, max {temp_in[var].max().values}, mean {temp_in[var].mean().values}")
    
    # for var in grid.data_vars:
    #     l.info(f"debug {var}: min {grid[var].min().values}, max {grid[var].max().values}, mean {grid[var].mean().values}")
    
    grid = grid[variables]
    
    l.info(f"Calculating PET with {args.method}")
    pet = calc_pet(
        grid,
        temp=temp_in,
        dem_model=DEM,
        method=args.method,
        press_correction=False,
        wind_correction=False,
        wind_altitude=False,
        freq='D',
        logger=l
        )
    
    l.info(f"PET calculated")
    
    year = pet.time.dt.year[0].item()
    
    fig, ax = plt.subplots()
    ax.set_title(f"Mean global penman-monteith PET for {year}")
    pet.mean(dim='time').plot(ax=ax, vmin=0, cmap='viridis')
    os.makedirs(f"data/4-output/PET_global/annual_plots", exist_ok=True)
    plt.savefig(f"data/4-output/PET_global/annual_plots/pet_{year}.png")
    plt.close()
    
    ds_out = xr.Dataset({'pet': pet}, coords=pet.coords, attrs=pet.attrs)

    if 'latitude' in ds_out.coords:
        ds_out = ds_out.rename({'latitude': 'y', 'longitude': 'x'})
    
    
    if args.clusterid is not None:  
        outdir = os.path.join(mod.root, "inmaps")
        os.makedirs(outdir, exist_ok=True)
        outdir = os.path.join(outdir, f'inmaps_pet_{args.tpf}_{args.method}_daily_*.nc')
    else:
        outdir = os.path.join("data/4-output/PET_global", f'{args.tpf}_{args.method}_daily_*.nc')
        os.makedirs(os.path.dirname(outdir), exist_ok=True)
    
    write_forcing(l,
                  ds_out, 
                  outdir, 
                  freq_out='Y', 
                  chunksize=1, 
                  decimals=2, 
                  time_units="days since 1900-01-01T00:00:00")

if __name__ == "__main__":
    l = setup_logging('data/0-log', '00_PET.log')
    l.info("Logging 00_PET.py to data/0-log/00_PET.log")
    try:
        parser = argparse.ArgumentParser(description='Build for PET globally, or for the region of a Wflow model.')
        parser.add_argument('--tpf', type=str, help='the temperature forcing', default='era5_daily_zarr')
        parser.add_argument('--method', type=str, help='the method to use for PET calculation', default='penman-monteith_tdew')
        parser.add_argument('--clusterid', type=str, help='the cluster id', default=None)
        parser.add_argument('--bbox', type=str, help='the bounding box for the region of interest', default=None)
        parser.add_argument('--tmin', type=str, help='the start time for the PET forcing %Y-%m-%d', default='1980')
        parser.add_argument('--tmax', type=str, help='the end time for the PET forcing %Y-%m-%d', default=None)
        args = parser.parse_args()

        main(l, args)
        
    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.info(traceback.format_exc())
        l.info("Exiting")
        exit(1)
