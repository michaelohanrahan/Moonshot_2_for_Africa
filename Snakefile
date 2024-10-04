import sys
import geopandas as gpd
import numpy as np
from pathlib import Path
import scripts
from scripts.helper import syscheck
from scripts.helper import create_directories
import os

#TODO: finalize the workflow and how its called
#      User defines a region or location --> builds a snake config --> runs build and instates for the region
#      ARGS: REGION, DATE ..... (¿ VIZ, PLOT, SAVE ?)

configfile: "config/snakeConfig.yml"

#return drive as either p: or /p/
DRIVE=syscheck()

# set the working directory
os.chdir(Path(f'{DRIVE}/moonshot2-casestudy/Wflow/africa'))

dirs = create_directories(config)
external, interdir, indir, outdir = dirs['external'], dirs['interim'], dirs['input'], dirs['output']

# since the wildcard is not known at the start of the pipeline, we need to define a function to get the clusters
def get_clusters(wildcards):return list(set(np.loadtxt(Path('data/2-interim/cluster_list.txt')).astype(int).tolist()))

rule all: 
    input: 
        expand(Path(outdir, 'models', '{cluster}', 'staticmaps', 'staticmaps.nc'), cluster=get_clusters(wildcards)),


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

rule clusterbasins:
    input: 
        basin_geojson = config['data']['basins']
    params: 
        crs = config['system']['crs'],
        method = config['methods']['cluster'],
        touches = config['methods']['touches'],
        plot = False,
        savefig = True,
        test_list = None, #can be None or list
        fill_rings = False
    output:models = Path(outdir, 'models', '{cluster}')#cluster=f'{{cluster}}')
    localrule: True
    script:
        "scripts/01_cluster_basins.py"

#build models
rule build_models:
    input:
        basin_geojson = 'data/2-interim/dissolved_basins.geojson',
        model = Path(outdir, 'models', '{cluster}')
    output:
        grid = Path(outdir, 'models','{cluster}','staticmaps', 'staticmaps.nc')
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

rule build_forcing:
    input:
        grid = Path(outdir, 'models','{cluster}','staticmaps', 'staticmaps.nc')
    params:
        tpf = config['data']['temp_precip_forcing'],
        method = config['methods']['pet_method'],
        tmin = config['methods']['tmin'], 
        tmax = config['methods']['tmax'],
    output:
        forcing = Path(outdir, 'models', '{cluster}', 'inmaps', 'inmaps_pet_prec_{params.tpf}_{params.method}_daily_*.nc')
    localrule: False
    shell:
        '''python scripts/00_forcing --tpf {params.tpf} --method {params.method} --tmin {params.tmin} --tmax {params.tmax}'''
    
rule run_wflow_forecast:
'''
:: Run a Wflow forecast for a given location and timestamp

args:
    cluster_selected: int, cluster that we want to run
returns:
    forcing.nc: netcdf file with the forcing data, written per year as '{model_root}/inmaps/inmaps_pet_prec_{tpf}_{.method}_daily_*.nc'
'''
    input:
        cluster_selected = 
    params:
    output:
    shell:

rule run_wflow_state:
    input:
    params:
    output:
    shell:

rule run_wflow_hist:
    input:
    params:
    output:
    shell: