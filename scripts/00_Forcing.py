import hydromt 
import xarray as xr 
import argparse
import traceback 
from scripts.helper import setup_logging
from shapely import geometry
import matplotlib.pyplot as plt
from hydromt_wflow import WflowModel as WM
from hydromt import DataCatalog
import os
import sys
from pathlib import Path
import pandas as pd
import geopandas as gpd


os.chdir(Path(r'p:\moonshot2-casestudy\Wflow\africa'))
l = setup_logging('data/0-log', '00_forcing.log')
l.info("building model assuming access to deltares_data catalog")
try:
    parser = argparse.ArgumentParser(description='Build Forcing for a Bounding Box')
    parser.add_argument('--tpf', type=str, help='the temperature forcing', default='era5_daily_zarr')
    parser.add_argument('--tpr', type=str, help='the precipitation forcing', default='era5_daily_zarr')
    parser.add_argument('--method', type=str, help='the method to use for PET calculation', default='penman-monteith_tdew')
    parser.add_argument('--clusterid', type=str, help='the cluster id', default='default_model')
    parser.add_argument('--tmin', type=str, help='the start time for the PET forcing %Y-%m-%d', default=None)
    parser.add_argument('--tmax', type=str, help='the end time for the PET forcing %Y-%m-%d', default=None)
    args = parser.parse_args()

    if args.tpf is None:
        l.error("Temperature forcing must be provided")
        raise ValueError("Temperature forcing must be provided")

    if args.method is None:
        l.error("Method must be provided")
        raise ValueError("Method must be provided")

    if args.clusterid is None:
        l.error("Cluster id must be provided")
        raise ValueError("Cluster id must be provided")

    root = f'src/3-model/wflow_build/{args.clusterid}'
    method = args.method
except SystemExit:
    l.error("Error parsing arguments")
    root = f'src/3-model/wflow_build/1844'
    method = 'penman-monteith_tdew'

config = f"wflow_sbm.toml"

mod = WM(root=root,
        mode='r',
        config_fn=config,
        data_libs=['p:/wflow_global/hydromt/deltares_data.yml'],
        logger=l)

bbox = mod.region.geometry
bbox_b = bbox.buffer(0.5).bounds.values[0]

dc=DataCatalog('p:/wflow_global/hydromt/deltares_data.yml')

xmin,ymin,xmax,ymax=bbox_b
if ymin > ymax:
    grid = dc.get_rasterdataset('era5_daily_zarr').sel(longitude=slice(xmin,xmax), latitude=slice(ymin,ymax)).chunk({"time":1})
    DEM = dc.get_rasterdataset('era5_orography').sel(longitude=slice(xmin,xmax), latitude=slice(ymin,ymax))
elif ymin < ymax:
    grid = dc.get_rasterdataset('era5_daily_zarr').sel(longitude=slice(xmin,xmax), latitude=slice(ymax,ymin)).chunk({"time":1})
    DEM = dc.get_rasterdataset('era5_orography').sel(longitude=slice(xmin,xmax), latitude=slice(ymax,ymin)).isel(time=0)

assert DEM.raster.identical_grid(grid) == True, "The grids are not identical"

t_present = pd.to_datetime(grid.time.values[-1])


variables = ["temp"]
if method == "debruin":
    variables += ["press_msl", "kin", "kout"]
elif method == "makkink":
    variables += ["press_msl", "kin"]
elif method == "penman-monteith_rh_simple":
    variables += ["temp_min", "temp_max", "rh", "kin"]
elif method == "penman-monteith_tdew":
    variables += [
        "temp_min",
        "temp_max",
        "wind10_u",
        "wind10_v",
        "temp_dew",
        "kin",
        "press_msl",
    ]

if 'u10' in grid.data_vars:
    grid = grid.rename({'u10':'wind10_u', 'v10':'wind10_v'})
if 'd2m' in grid.data_vars:
    grid = grid.rename({'d2m':'temp_dew'})

P = grid['precip']
precip = hydromt.workflows.forcing.precip(P,da_like=P, freq='D', logger=l).sel(time=slice('1980',t_present))

# grid = grid.rename({'u10':'wind_u', 'v10':'wind_v'})
temp_in = hydromt.workflows.forcing.temp(
    grid['temp'], 
    dem_model=DEM, 
    dem_forcing=DEM, 
    lapse_correction=False, 
    logger=l, 
    freq=None
    )

if (
    "penman-monteith" in method
):
    l.info("Calculating max and min temperature")
    temp_max_in = hydromt.workflows.forcing.temp(
        grid['temp_max'],
        dem_model=DEM,
        dem_forcing=DEM,
        lapse_correction=False,
        logger=l,
        freq=None
    )
    temp_max_in = temp_max_in.rename('temp_max')
    temp_min_in = hydromt.workflows.forcing.temp(
        grid['temp_min'],
        dem_model=DEM,
        dem_forcing=DEM,
        lapse_correction=False,
        logger=l,
        freq=None
    )
    temp_min_in = temp_min_in.rename('temp_min')
    temp_in=xr.merge([temp_in, temp_max_in, temp_min_in])
    
else:
    temp_in = hydromt.workflows.forcing.temp(T, dem_model=DEM, dem_forcing=DEM, lapse_correction=False)

pet = hydromt.workflows.forcing.pet(
    grid[variables[1:]], 
    temp=temp_in, 
    dem_model=DEM, 
    method=method, 
    press_correction=False,
    wind_correction=False,
    wind_altitude=False,
    freq='D',
    logger=l)

if "penman-monteith" in method:
    temp_in = temp_in.drop(['temp_max', 'temp_min'])

temp_out=hydromt.workflows.forcing.resample_time(
    temp_in, 
    freq='D', 
    upsampling="bfill",  # we assume right labeled original data
    downsampling="mean",
    label="right",
    closed="right",
    conserve_mass=False,
    logger=l,
)

mod.forcing['pet'] = pet
mod.forcing['temp'] = temp
mod.forcing['precip'] = precip
os.makedirs(os.path.join(mod.root, "inmaps"), exist_ok=True)
mod.config['input']['path_forcing'] = f'inmaps/inmaps_pet_prec_{tpf}_{method}_daily*.nc'
mod.write_config()
mod.write_forcing(freq_out='Y')