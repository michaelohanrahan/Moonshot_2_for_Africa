import os
import geopandas as gpd

os.chdir(r'p:\moonshot2-casestudy\Wflow\africa')

rule clusterbasins:
    input: 
        basin_geojson = r"data\2-interim\GIS\basins_mainland_and_madagascar.geojson"
    params: 
        method = "domain_method",
        intersect_method = "centroid",
        crs = "EPSG:4326",
        plot = False,
        savefig = False,
        test_list = None, #can be None or list
        fill_rings = False
    output:
        cluster_out = os.path.join('data', '2-interim', 'dissolved_basins.geojson')
    script:
        "scripts/01_cluster_basins.py"

#Next rule uses cluster list
rule derive_basins:
    input:


