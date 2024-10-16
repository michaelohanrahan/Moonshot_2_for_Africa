import os
import sys
import geopandas as gpd
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from pathlib import Path
import argparse
from icecream import ic

def main(root:str, 
         gauges, 
         new_root:str = None, 
         mode:str = "w", 
         basename:str = "hourly", 
         index_col:str = "id",
         snap_to_river:bool = True,
         max_dist:int = 1000,
         derive_subcatch:bool = True,
         crs='epsg:4326',
         config_old="wflow_sbm_template.toml",
         config_new="wflow_sbm_add_gauges.toml",
         ignore_list=None,
         **kwargs):
    
    """
    SUMMARY:
    """
    if sys.platform == "win32":
        DRIVE='p:'
    else:
        DRIVE='/p'
    
    if not Path(gauges).is_absolute():
        gauges = Path(Path.cwd(), gauges)
        
    gauges = gpd.read_file(gauges, 
                           crs=crs, 
                           read_geometry=False, 
                           columns=[index_col, 'x', 'y'])
    # ic(gauges)
    if ignore_list:
        gauges = gauges[~gauges[index_col].isin(ignore_list)]
    # ic(gauges)
    
    if 'geometry' not in gauges.columns:
        if 'x' in gauges.columns and 'y' in gauges.columns:
            gauges['geometry'] = gpd.points_from_xy(gauges['x'], gauges['y'])
        elif 'lon' in gauges.columns and 'lat' in gauges.columns:
            gauges['geometry'] = gpd.points_from_xy(gauges['lon'], gauges['lat'])
        else:
            raise ValueError('Gauges file does not contain geometry, x/y or lon/lat columns')
        
    gauges = gauges.set_geometry(col='geometry')
    gauges = gauges.set_crs(crs)  # Explicitly set the CRS for the GeoDataFrame
    
    if new_root is None:
        new_root = root
        mode = "w+"

    logger = setuplog("build", log_level=20)
    #ic(root)
    #ic(config_old)
    w = WflowModel(
        root=root,
        mode="r",
        config_fn = config_old,
        data_libs = [f'{DRIVE}/wflow_global/hydromt/deltares_data.yml'],
        logger=logger,
        )
    logger.info(f"Reading model from {root}")
    w.read()
    
    logger.info(f"Setting CRS to {crs}")
    w.set_crs(crs)
    w.set_root(
        root=new_root,
        mode=mode
        )
    
    w.set_crs(crs)
    # ic(w.crs)
    w.setup_gauges(
        gauges_fn=gauges,
        snap_to_river=snap_to_river,
        derive_subcatch=derive_subcatch,
        index_col=index_col,
        basename=basename,
        max_dist=max_dist,
        toml_output='csv',
        gauge_toml_header = ["Q", "P", "PET"],
        gauge_toml_param = ["lateral.river.q_av", "vertical.precipitation", "vertical.potential_evaporation"],
        **kwargs
    )
    
    
    w.config['input']['path_static'] = f'staticmaps/staticmaps.nc'
    print('writing config')
    w.write_config(config_name=config_new)
    print('writing geoms')
    w.write_geoms()
    print('writing grid')
    w.write_grid()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update wflow gauges configuration.")
    parser.add_argument("cwd", help="Current working directory", type=str)
    parser.add_argument("config_root", help="Root directory for configuration", type=str)
    parser.add_argument("gauges", help="Gauges file", type=str)
    parser.add_argument("--ignore_list", help="string to return separated list of integers", type=str, default="data/2-interim/ignore_list.txt") 
    parser.add_argument("--new_root", help="New root directory (optional)", type=str, default=None)
    parser.add_argument("--mode", help="Mode for opening the model (r, w, w+)", type=str, default="w")
    parser.add_argument("--basename", help="Basename for the gauges", type=str, default="hourly")
    parser.add_argument("--index_col", help="Index column for the gauges", type=str, default="id")
    parser.add_argument("--snap_to_river", help="Snap gauges to river", type=bool, default=True)
    parser.add_argument("--max_dist", help="Maximum distance for snapping", type=int, default=1000)
    parser.add_argument("--derive_subcatch", help="Derive subcatchments", type=bool, default=True)
    parser.add_argument("--crs", help="Coordinate reference system", type=str, default="epsg:4326")
    parser.add_argument("--config_old", help="Configuration file", type=str, default="wflow_sbm_template.toml")
    parser.add_argument("--config_new", help="Configuration file", type=str, default="wflow_sbm_add_gauges.toml")
    parser.add_argument("--kwargs", help="Additional keyword arguments", type=dict, default={})


    args = parser.parse_args()
    os.chdir(args.cwd)
    root = os.getcwd()
    
    # #ic(root)
    
    if args.new_root:
        new_root = os.path.join(root, args.new_root)
    else:
        raise ValueError("New root directory must be provided.")
    
    if args.ignore_list and os.path.exists(args.ignore_list):
        
        with open(args.ignore_list, 'r') as file:
            ignore_list = [int(line.strip()) for line in file if line.strip().isdigit()]
        #ic(ignore_list)
    else:
        ignore_list = None
        
    root = os.path.join(root, args.config_root)
    # #ic(ignore_list)
    main(root=root, 
         gauges=args.gauges, 
         new_root=new_root,
         mode=args.mode,
         basename=args.basename,
        index_col=args.index_col,
        snap_to_river=args.snap_to_river,
        max_dist=args.max_dist,
        derive_subcatch=args.derive_subcatch,
        crs=args.crs,
        config_old=args.config_old,
        config_new=args.config_new,
        ignore_list=ignore_list,
        kwargs=args.kwargs
    )