import os

cluster_method = "domain_method" #{domain_method, ...}
crs = "EPSG:4326"
basin_path = r"data\2-interim\GIS\basins_mainland_and_madagascar.geojson"

rule all:
    input: 
        cluster_out = os.path.join('data', '2-interim', 'clustered', cluster_method+'.geojson')

rule cluster:
	input:
		basin_geojson = basin_path
	params:
		method = cluster_method,
		crs = crs
	script:
		"scripts/01_cluster_basins.py"
	output:
		out = os.path.join('data', '2-interim', 'clustered', cluster_method+'.geojson')