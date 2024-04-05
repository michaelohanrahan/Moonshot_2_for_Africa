import os

cluster_method = "domain_method"
crs = "EPSG:4326"
basin = r"p:\moonshot2-casestudy\Wflow\africa\data\2-interim\GIS\basins_mainland_and_madagascar.geojson"

rule all:
	os.path.join('data', '2-interim', 'clustered', cluster_method+'.geojson')

rule cluster:
	input:
		basin_geojson = basins
	params:
		method = cluster_method
		crs = crs
	output:
		out = os.path.join('data', '2-interim', 'clustered', cluster_method+'.geojson')