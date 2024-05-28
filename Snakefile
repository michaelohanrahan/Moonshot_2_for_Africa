import sys
import geopandas as gpd
import numpy as np
from pathlib import Path
import scripts
from scripts.helper import syscheck
from scripts.helper import create_directories
import os

configfile: "snakeConfig.yml"

DRIVE=syscheck()

# set the working directory
os.chdir(Path(f'{DRIVE}/moonshot2-casestudy/Wflow/africa'))

#Setup the data directories
external = Path('data/1-external')
external.mkdir(exist_ok=True, parents=True)
#This directory will require some manual intervention: 
# !looking here for a basin file to cluster

interdir = Path('data/2-interim')
interdir.mkdir(exist_ok=True, parents=True)

indir = Path('data/3-input')
indir.mkdir(exist_ok=True, parents=True)

outdir = Path('data/4-output')
outdir.mkdir(exist_ok=True, parents=True)

dirs = create_directories(config)
external, interdir, indir, outdir = dirs


rule all:
    input: expand(Path(outdir, model, '{cluster}', staticmaps, 'staticmaps.nc'), cluster=clusters)

# Cluster basins
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
    output:
        cluster_out = Path('data/2-interim/dissolved_basins.geojson'),
        clusters = np.load_txt(Path('data/2-interim/cluster_list.txt')).astype(int).tolist().unique()
    script:
        "scripts/01_cluster_basins.py"

# TODO: add the building rule
rule build_models:
    input:
        basin_geojson = 'data/2-interim/dissolved_basins.geojson'
    output:
         expand(Path(outdir, model, '{cluster}', staticmaps, 'staticmaps.nc'), cluster=clusters)
    script:
        "scripts/02_hydromt_wflow_build.py"