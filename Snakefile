import os

cluster_method = "domain_method" #{domain_method, ...}
crs = "EPSG:4326"
basin_path = r"data\2-interim\GIS\basins_mainland_and_madagascar.geojson"

rule all:
    input: 
        cluster_out = os.path.join('data', '2-interim', 'clustered', cluster_method+'.geojson')

rule clusterbasins:
    input: 
        basin_path
    params: 
        method = cluster_method,
        crs = crs
    output:
        os.path.join('data', '2-interim', 'clustered', cluster_method+'.geojson')
	script:
		"scripts/01_cluster_basins.py"