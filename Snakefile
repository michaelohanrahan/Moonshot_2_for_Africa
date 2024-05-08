import geopandas as gpd
from pathlib import Path
from scripts.helper import syscheck

DRIVE=syscheck()

# set the working directory
os.chdir(Path(f'{DRIVE}/moonshot2-casestudy/Wflow/africa'))

rule all:
    input:
        'data/2-interim/dissolved_basins.geojson'

# Cluster basins
#TODO: fix endhoreic basins
rule clusterbasins:
    input: 
        basin_geojson = config['files']['basins']
    params: 
        crs = config['system']['crs'],
        method = config['methods']['cluster'],
        touches = config['methods']['touches'],
        plot = False,
        savefig = False,
        test_list = None, #can be None or list
        fill_rings = False
    output:
        cluster_out = Path('data/2-interim/dissolved_basins.geojson')
    script:
        "scripts/01_cluster_basins.py"
