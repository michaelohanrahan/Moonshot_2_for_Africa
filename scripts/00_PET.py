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


def init_WM(l,
            root: str,
            config: str)->WM:
    '''
    initialize a new wflow model
    '''
    mod = WM(root=root, 
             mode='w+',
             config_fn=config,
             data_libs='deltares_data',
             logger=l)
    
    return mod

def build_WM(l, 
             mod):
    '''
    '''
    

def main(l,
         tpf: str,
         bbox: list | tuple,
         method: str, 
         mod: WM,
         clusterid: int,
         Test: bool = None)-> None:
    '''
    Takes a bounding box and extracts
    '''
    mod.read_config()
    mod.read_grid()
    
    #:: Assuming full time bounds 
    dc = DataCatalog('deltares_data')
    ds = dc.get_rasterdataset(tpf)
    tmin = pd.Timestamp(ds.time.min().values).strftime('%Y-%m-%dT%H:%M:%S')
    tmax = pd.Timestamp(ds.time.max().values).strftime('%Y-%m-%dT%H:%M:%S')
    
    #:: LOGGING
    l.info(f"creating the PET forcing for model {mod.root}")
    l.info(f"using the {method} method")
    l.info(f"using data source {tpf}")
    
    #:: CONFIG for time
    mod.config['calendar'] = 'proleptic_gregorian'
    l.info(f"calendar set to {mod.config['calendar']}")
    mod.config['starttime'] = tmin
    l.info(f"start time set to {mod.config['starttime']}")
    mod.config['endtime'] = tmax
    l.info(f"end time set to {mod.config['endtime']}")
    mod.config['time_units'] = 'days since 1900-01-01 00:00:00'
    l.info(f"time units set to {mod.config['time_units']}")
    mod.config['timestepspecs'] = 86400
    l.info(f"timestep set to {mod.config['timestepspecs']}")
    mod.config['log'] = os.path.relpath(f'data/0-log/{clusterid}/PET.log', mod.root)
    l.info(f"log file set to {mod.config['log']}")
    mod.config['input']['path_forcing'] = f'inmaps_pet_{tpf}_{method}_daily.nc'
    
    mod.write_config()
    
    l.info(f"Setting up the PET forcing for model {mod.root}")
    mod.setup_temp_pet_forcing(temp_pet_fn=tpf,
                               pet_method=method,
                               press_correction=True,
                               temp_correction=True,
                               wind_correction=True,
                               wind_altitude=10, #default
                               chunksize=1
                               )
    l.info(f"finished PET setup")
    l.info(f"writing the PET forcing to {mod.config['input']['path_forcing']}")
    
    mod.write_forcing(freq_out='Y')
        
    
if __name__ == "__main__": 
    os.chdir(Path(r'p:\moonshot2-casestudy\Wflow\africa'))
    l = setup_logging('data/0-log', '00_PET.log')
    l.info("building model assuming access to deltares_data catalog")
    try:
        if 'ipykernel' in sys.modules:
            l.info("Running in Jupyter")
            
            bbox = [34.33,-20.12,34.95,-19.30] #xmin, ymin, xmax, ymax
            root = 'data/2-interim/'
            method = 'penman-monteith_tdew'
            bm_toml = 'base_model.toml'
            bm_root = 'data/2-interim'
            clusterid = 0000
            root = Path(root, str(clusterid)).as_posix()
            config = f"wflow_sbm_{method}_{clusterid}.toml"
            tpf = 'era5_daily_zarr'
            
            if not os.path.exists(Path(bm_root,bm_toml)):
                mod = init_WM(l,root,config)
                mod.setup_basemaps(region={'bbox': list(bbox)}, 
                                   hydrography_fn='merit_hydro',
                                   res=0.008333333333333333,)
                mod.write_grid()
                mod.config['case'] = f'wflow_sbm model for creating PET forcing for East Africa cluster {clusterid}'
                mod.write_config()
            else:
                mod = WM(root=root,
                            mode='w+',
                            config_fn=config,
                            logger=l)
            #setup PET for the model
            main(l,tpf,bbox,method,mod,clusterid,Test=True)
        else:     
            parser = argparse.ArgumentParser(description='Process some integers.')
            parser.add_argument('bbox', metavar='N', type=float, nargs='+', help='an integer for the accumulator')
            args = parser.parse_args()
            
            if len(args.bbox) != 4 or args.bbox is None:
                l.error("Bounding box must have 4 values in a list or tuple")
                raise ValueError("Bounding box must have 4 values in a list or tuple")
            l.info(f"Bounding box: {args.bbox}")
            main(l,list(args.bbox))
    
    except Exception as e:
        l.error(f"Error: {e}")
        l.error(traceback.format_exc())
        sys.exit()