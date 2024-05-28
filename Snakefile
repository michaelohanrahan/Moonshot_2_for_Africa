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