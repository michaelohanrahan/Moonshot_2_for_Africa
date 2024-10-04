import sys
import geopandas as gpd
import numpy as np
from pathlib import Path
import scripts
from helper import syscheck
from helper import create_directories
import os

#TODO: finalize the workflow and how its called
#      User defines a region or location --> builds a snake config --> runs build and instates for the region
#      ARGS: REGION, DATE ..... (¿ VIZ, PLOT, SAVE ?)

configfile: "config/snakeConfig.yml"

#return drive as either p: or /p/
DRIVE=syscheck()

# set the working directory
os.chdir(Path(f'{DRIVE}/moonshot2-casestudy/Wflow/africa'))

data_dirs = create_directories(config)
external, interdir, indir, outdir = data_dirs['external'], data_dirs['interim'], data_dirs['input'], data_dirs['output']

source_wflow = config['data']['models']['wflow']
source_sfincs = config['data']['models']['sfincs_boundary']

clusters = np.loadtxt(Path('data/2-interim/cluster_list_1844.txt')).astype(int).tolist()

rule all: 
    input: 
        expand(Path(source_wflow, '{cluster}', 'staticmaps.nc'), cluster=clusters),
        expand(Path(source_sfincs, '{cluster}', 'staticmaps.nc'), cluster=clusters),
        expand(Path(source_sfincs, '{cluster}', 'wflow_sbm.toml'), cluster=clusters)

'''
:: Cluster the basins based on the given method

args:
    basin_geojson: str, path to the basin GeoJSON file
    crs: str, crs of the basin GeoJSON file
    method: str, method to use for clustering
    touches: bool, if True, the basins will be clustered based on the touching basins
    plot: bool, if True, the clusters will be plotted
    savefig: bool, if True, the plot will be saved
    test_list: list, list of basins to test the clustering on
    fill_rings: bool, if True, the rings of the basins will be filled

#TODO:  Think about how the clustering actually fits into the workflow.. 
        Since this is the first step below, but the indexes of the clustering are already known for the wildcards. 
'''

# rule clusterbasins:
#     input: 
#         basin_geojson = config['data']['basins']
#     params: 
#         crs = config['system']['crs'],
#         method = config['methods']['cluster'],
#         touches = config['methods']['touches'],
#         plot = False,
#         savefig = True,
#         test_list = None, #can be None or list
#         fill_rings = False
#     output:models = Path(outdir, 'models', '{cluster}')#cluster=f'{{cluster}}')
#     localrule: True
#     script:
#         "scripts/01_cluster_basins.py"

#build models
rule build_models:
    input:
        basin_geojson = 'data/2-interim/dissolved_basins.geojson',
        model = Path(source_wflow, '{cluster}', "staticmaps.nc")
    output:
        grid = Path(source_wflow, 'models','{cluster}','staticmaps', 'staticmaps.nc')
    localrule: False #Determine the necessary recources by RAM demand and add to the profile
    script:
        "scripts/02_hydromt_wflow_build.py"

'''
:: Build the forcing file for every created model

The forcing file is lightweight, containing temperature, precipitation, and PET data for the model. 
The data is written per year as "forcing_{year}.nc".

args:
    tpf: str, name of the temp_precip_forcing datasource in the datacatalog
    method: str, method to use for calculating PET
    tmin: str, year to start from
    tmax: str, year to end at (None defaults to most recent, present year)
returns:
    forcing.nc: netcdf file with the forcing data, written per year as '{model_root}/inmaps/inmaps_pet_prec_{tpf}_{.method}_daily_*.nc'

#TODO: IDEALLY we no longer build forcing, we run from DD catalog checking for instates and data availability 
        including download and processing of the data.

'''

# rule build_forcing:
#     input:
#         grid = Path(outdir, 'models','{cluster}','staticmaps', 'staticmaps.nc')
#     params:
#         tpf = config['data']['temp_precip_forcing'],
#         method = config['methods']['pet_method'],
#         tmin = config['methods']['tmin'], 
#         tmax = config['methods']['tmax'],
#     output:
#         forcing = Path(outdir, 'models', '{cluster}', 'inmaps', 'inmaps_pet_prec_{params.tpf}_{params.method}_daily_*.nc')
#     localrule: False
#     shell:
#         '''python scripts/00_forcing --tpf {params.tpf} --method {params.method} --tmin {params.tmin} --tmax {params.tmax}'''


'''
create the gaugemap from your preferred geojson file of gauges... best to make this with some manual checks to ensure the right gauges are being used
1: src/pre/assess_discharge_data.py to perform health checks and gather all gauges to be used (assuming all data is collected in a datacatalog)
2: scripts\pixi_run_gaugemap.bat will build a gaugemap from the output geojson and the dependency graph
 -- you can recursively modify the ignorelist in the interim folder to remove gauges to simplify the graph
    #TODO: automatically optimise the graph
3: src/pre/create_discharge_data.py to to create the discharge dataset from the gauges geojson
adding default target to false to encourage snakemake to run these rules before defining the wildcards
'''


rule create_gaugemap:
    input:
        gauges = Path(interdir, "model_inflow_points", '{cluster}', "gauges.geojson")
    output:
        gridfile = Path(source_sfincs, '{cluster}', "staticmaps.nc"),
        toml = Path(source_sfincs, '{cluster}', "wflow_sbm.toml")
    params:
        cwd=Path(os.getcwd()).as_posix(),
        config_root=Path(source_wflow, '{cluster}'),
        new_root=Path(source_sfincs, '{cluster}'),
        mode="w+",
        ignorelist=Path(interdir, "model_inflow_points", '{cluster}', "ignore_list.txt"),
        basename="sfincs",
        index_col="wflow_id",
        max_dist=1000,
        crs=config['system']['crs'],
        config_old="wflow_sbm.toml",
        config_new="wflow_sbm.toml"
    localrule: True
    cores: 1
    shell:
        """
        pixi run python src/pre/update_gauges.py {params.cwd} {params.config_root} {input.gauges} \
            --new_root "{params.new_root}" \
            --mode "{params.mode}" \
            --basename "{params.basename}" \
            --index_col "{params.index_col}" \
            --ignore_list "{params.ignorelist}" \
            --snap_to_river True \
            --max_dist {params.max_dist} \
            --derive_subcatch True \
            --crs "{params.crs}" \
            --config_old "{params.config_old}" \
            --config_new "{params.config_new}"
        """