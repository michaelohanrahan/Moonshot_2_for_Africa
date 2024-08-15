import sys
import geopandas as gpd
import numpy as np
from pathlib import Path
import scripts
from scripts.helper import syscheck
from scripts.helper import create_directories
import os

configfile: "config/snakeConfig.yml"

DRIVE=syscheck()

# set the working directory
os.chdir(Path(f'{DRIVE}/moonshot2-casestudy/Wflow/africa'))

dirs = create_directories(config)
external, interdir, indir, outdir = dirs['external'], dirs['interim'], dirs['input'], dirs['output']

# since the wildcard is not known at the start of the pipeline, we need to define a function to get the clusters
def get_clusters(wildcards):return list(set(np.loadtxt(Path('data/2-interim/cluster_list.txt')).astype(int).tolist()))

rule all: 
    input: 
        Path(outdir, 'models', '1716', 'staticmaps', 'staticmaps.nc')

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
    output:models = expand(Path(outdir, 'models', '{cluster}'), cluster=f'{{cluster}}')
    script:
        "scripts/01_cluster_basins.py"

#build models
rule build_models:
    input:
        basin_geojson = 'data/2-interim/dissolved_basins.geojson',
        model = lambda wildcards: expand(str(Path(outdir, 'models', '{cluster}')), cluster=get_clusters(wildcards))
    output:
        grid = Path(outdir, 'models','{cluster}','staticmaps', 'staticmaps.nc')
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
'''

rule build_forcing:
    input:
        grid = Path(outdir, 'models','{cluster}','staticmaps', 'staticmaps.nc')
    params:
        tpf = config['data']['temp_precip_forcing'],
        method = config['methods']['pet_method'],
        tmin = config['methods']['tmin'], #can be year
        tmax = config['methods']['tmax'], #can be None and defaults to most recent
    output:
        forcing = Path(outdir, 'models', '{cluster}', 'inmaps', 'inmaps_pet_prec_{params.tpf}_{params.method}_daily_*.nc')
    shell:
        '''python scripts/00_forcing --tpf {params.tpf} --method {params.method} --tmin {params.tmin} --tmax {params.tmax}'''