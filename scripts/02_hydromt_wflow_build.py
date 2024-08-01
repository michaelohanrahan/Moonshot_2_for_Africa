import datetime
import os
from hydromt_wflow import WflowModel
import hydromt
import geopandas as gpd

# setup logging
from hydromt.log import setuplog
logger = setuplog("Moonshot 2 - Africa", log_level=10)

# hard-coded input for testing
from sys import platform
if platform == "linux" or platform == "linux2":
    ROOT = "/p/moonshot2-casestudy/Wflow/africa/"
elif platform == "win32":
    ROOT = "p:/moonshot2-casestudy/Wflow/africa/"

INDEX_COL = "cluster_key"
CLUSTERED_GEOMETRIES = os.path.join(ROOT, "data/2-interim/dissolved_basins.geojson")

# global settings for Wflow model
MODE = "w"
BUILD_CONFIG = os.path.join(ROOT, "config/02_hydromt-build-full.yml")
FORCING_CONFIG = os.path.join(ROOT, "config/02_hydromt-update-era5.yml")
WFLOW_ROOT = os.path.join(ROOT, "src/3-model/wflow_build") #TODO model locations?

# snakemake input
# CLUSTERED_GEOMETRIES = snakemake.params.clustered_geometries

def log_basin_id(basin_id: int, log_file: str):
    with open(log_file, 'a') as file:
        log_entry = f"Failed basin id {basin_id}: {datetime.datetime.now()}\n"
        file.write(log_entry)

def get_catalog() -> str:
    from sys import platform
    if platform == "linux" or platform == "linux2":
        path_catalog = "/p/moonshot2-casestudy/Wflow/africa/config/deltares_data_custom.yml"
        # raise NotImplementedError("Linux data catalog not yet implemented in this script") TODO remove from code when it works properly
    elif platform == "win32":
        path_catalog = "p:/wflow_global/hydromt/deltares_data.yml"
    else:
        raise SystemError(f"operating system {platform} not supported by this script")
    return path_catalog

def get_root(index: int, prefix=None) -> str:
    root = f"{prefix}{index}" if prefix else f"{index}"
    full_root = os.path.join(WFLOW_ROOT, root)
    return full_root

def get_geom(index: int, all_geoms: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    geom = all_geoms[all_geoms[INDEX_COL] == index]
    geom = geom.explode() #MULTIPOLYGON to POLYGON, TODO remove when no longer necessary
    return geom

def create_model(root: str, geom: gpd.GeoDataFrame, uparea: float = 1, res: float = None) -> None:
    # get data catalog, either Linux or Windows
    data_libs = [get_catalog()]
    config = hydromt.config.configread(config_fn=BUILD_CONFIG)
    # update model resolution if provided
    if res:
        config['setup_basemaps']['res'] = res
    # initialize model in root directory linked to basin ID
    w = WflowModel(root=root, mode=MODE, data_libs=data_libs, logger=logger)
    # see for region settings: https://deltares.github.io/hydromt/latest/user_guide/model_region.html
    region = {"basin": geom, "uparea": uparea}
    # build model using geometry file linked to basin ID
    w.build(region=region, write=True, opt=config)
    return w

if __name__ == "__main__":
    log_failed = os.path.join(WFLOW_ROOT, "failed_builds.log")
    with open(log_failed, 'w') as file:
        file.write('')

    logger.info(f"CUSTOM: reading geometry file with all basins")
    clustered_basins = gpd.read_file(CLUSTERED_GEOMETRIES)
    basin_ids = list(clustered_basins[INDEX_COL])
    logger.info(f"CUSTOM: geometry file contains {len(basin_ids)} basins")
    for i, basin_id in enumerate(basin_ids):
        logger.info(f"CUSTOM: processing basin id {basin_id} ({i + 1}/{len(basin_ids)})")
        logger.info(f"CUSTOM: creating root linked to basin id {basin_id}")
        root = get_root(index=basin_id)
        logger.info(f"CUSTOM: reading geometry file and finding geom with index {basin_id}")
        geom = get_geom(index=basin_id, all_geoms=clustered_basins)
        try:
            logger.info(f"CUSTOM: start initializing and building Wflow model with index {basin_id}")
            create_model(root=root, geom=geom)
        except Exception as e:
            logger.error(f"Error processing basin id {basin_id}: {e}")
            log_basin_id(basin_id)
