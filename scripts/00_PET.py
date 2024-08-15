import hydromt 
import xarray as xr 
import argparse
import traceback 
from helper import setup_logging
from shapely import geometry
import matplotlib.pyplot as plt
from hydromt_wflow import WflowModel as WM
from hydromt import DataCatalog
import os
import sys
from pathlib import Path
import pandas as pd
import geopandas as gpd


def init_WM(l,
            root: str,
            config: str)->WM:
    '''
    initialize a new wflow model
    '''
    mod = WM(root=root, 
             mode='w+',
             config_fn=config,
             data_libs=['p:/wflow_global/hydromt/deltares_data.yml'],
             logger=l)
    
    return mod

    

def main(l,
         root: str,
         tpf: str,
         method: str, 
         mod: WM,
         clusterid: int,
         tmin:str=None,
         tmax:str=None)-> None:
    '''
    Takes a bounding box and extracts
    '''
    mod.read_config()
    mod.read_grid()
    mod.read_geoms()
    
    #:: Assuming full time bounds 
    dc = DataCatalog('deltares_data')
    l.info([key for key in dc.keys if 'era5' in key])
    l.info(dc.get_rasterdataset(tpf))
    ds = dc.get_rasterdataset(tpf)
    
    if tmin is None:
        tmin = pd.Timestamp(ds.time.min().values).strftime('%Y-%m-%dT%H:%M:%S')
    else:
        tmin = pd.Timestamp(tmin).strftime('%Y-%m-%dT%H:%M:%S')
    if tmax is None:
        tmax = pd.Timestamp(ds.time.max().values).strftime('%Y-%m-%dT%H:%M:%S')
    else:
        tmax = pd.Timestamp(tmax).strftime('%Y-%m-%dT%H:%M:%S')
        
    l.info(f"Time bounds: {tmin} to {tmax}")
    
    #:: LOGGING
    l.info(f"creating the PET forcing for model {mod.root}")
    l.info(f"using the {method} method")
    l.info(f"using data source {tpf}")
    #:: CONFIG for time
    # mod.config['calendar'] = 'proleptic_gregorian'
    mod.config['starttime'] = tmin    
    mod.config['endtime'] = tmax
    mod.config['time_units'] = 'days since 1900-01-01 00:00:00'
    # mod.config['timestepspecs'] = 86400
    mod.config['log'] = os.path.relpath(f'data/0-log/{clusterid}/PET.log', mod.root)
    l.info(f"log file set to {mod.config['log']}")
    os.makedirs(os.path.join(mod.root, "PET"), exist_ok=True)
    forcing_path = os.path.relpath(mod.root, os.path.join(root, f'inmaps_pet_prec_{tpf}_{method}_daily*.nc'))
    mod.config['input']['path_forcing'] = forcing_path
    mod.write_config()
    
    l.info(f"Setting up the PET forcing for model {mod.root}")
    
    mod.setup_temp_pet_forcing(temp_pet_fn=ds,
                               pet_method=method,
                               press_correction=True,
                               temp_correction=True,
                               wind_correction=True,
                               wind_altitude=10, #default
                               chunksize=1
                               )
    
    l.info(f"finished PET setup")
    mod.setup_precip_forcing(precip_fn=ds, 
                             precip_clim_fn=None,
                             chunksize=1)
    
    l.info(f"forcing attrs: {mod.forcing.keys()}")
    
    mod.forcing['pet'].attrs['pet_fn'] = f"{tpf}_{method}"
    mod.forcing['temp'].attrs['temp_fn'] = f"{method}_corrected"
    mod.forcing['precip'].attrs['precip_fn'] = f"{tpf}"
    
    l.info(f"writing the forcing to {mod.root}")
    
    mod.write_forcing(freq_out='Y', 
                      )
    mod.write_geoms()
        
    
if __name__ == "__main__": 
    os.chdir(Path(r'p:\moonshot2-casestudy\Wflow\africa'))
    l = setup_logging('data/0-log', '00_PET.log')
    l.info("building model assuming access to deltares_data catalog")
    try:
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument('--tpf', type=str, help='the temperature forcing', default='era5_daily_zarr')
        parser.add_argument('--tpr', type=str, help='the precipitation forcing', default='era5_daily_zarr')
        parser.add_argument('--method', type=str, help='the method to use for PET calculation', default='penman-monteith_tdew')
        parser.add_argument('--clusterid', type=int, help='the cluster id', default='default_model')
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
        config = f"wflow_sbm.toml"

        mod = WM(root=root,
                    mode='w+',
                    config_fn=config,
                    data_libs=['p:/wflow_global/hydromt/deltares_data.yml'],
                    logger=l)
        
        mod.read_geoms()
        bbox = mod.geoms['basins'].total_bounds
        
        l.info(f"Model bounding box is {bbox}")
        tmp_root  = os.path.join(mod.root, 'tmp')
        os.makedirs(tmp_root, exist_ok=True)
        # if not os.path.exists(Path(bm_root,bm_toml)):
        
        mod = init_WM(l,tmp_root,config)
        
        mod.setup_basemaps(region={'bbox': list(bbox)}, 
                            hydrography_fn='merit_hydro',
                            res=0.25,)
        mod.write_grid()
        mod.config['case'] = f'wflow_sbm model for creating PET forcing for East Africa cluster {args.clusterid}'
        mod.write_config()
        
        main(l,
             root,
            args.tpf,
            args.method,
            mod,
            args.clusterid,
            args.tmin,
            args.tmax)
        
        
    
    except Exception as e:
        l.error(f"Error: {e}")
        l.error(traceback.format_exc())
        sys.exit()