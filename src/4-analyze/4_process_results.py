import os
import glob

import xarray as xr
import geopandas as gpd
import pandas as pd

from shapely.geometry.point import Point
from shapely.ops import nearest_points

path_forecast = "p:/moonshot2-casestudy/Wflow/workshop-prep/wflow_forecast/workshop_idai_chirps"
points_of_interest = "p:/moonshot2-casestudy/SFINCS/models/source_points_sfincs.geojson"

files = glob.glob(os.path.join(path_forecast, "**/output.nc"))
points = gpd.read_file(points_of_interest)

clusters = [1814, 1816, 1844]

SEARCH_RADIUS = 0.05 # degrees
VAR = "q_river"

df = pd.DataFrame()
for cluster in clusters:
    path_output = f"p:/moonshot2-casestudy/Wflow/workshop-prep/wflow_forecast/workshop_idai_chirps/{cluster}/output.nc"
    ds = xr.load_dataset(path_output)[VAR]

    path_wflow_river = f"p:/moonshot2-casestudy/Wflow/africa/data/3-input/wflow_build/{cluster}/staticgeoms/rivers.geojson"
    path_wflow_cluster = f"p:/moonshot2-casestudy/Wflow/africa/data/3-input/wflow_build/{cluster}/staticgeoms/basins.geojson"

    gdf_riv = gpd.read_file(path_wflow_river)
    gdf_cluster = gpd.read_file(path_wflow_cluster)
    
    points_cluster = gpd.sjoin(points, gdf_cluster, op='within')
    print(f"found {len(points_cluster)}/{len(points)} points of interest within cluster {cluster}")
    
    points_snapped = []
    for id, point in points_cluster.iterrows():

        circle = point.buffer(SEARCH_RADIUS)
        lines_within_circle = gdf_riv[gdf_riv.intersects(circle)]
            
        line_snapped = lines_within_circle.loc[lines_within_circle['strord'].idxmax()]
        point_snapped = nearest_points(point.geometry, line_snapped.geometry)[1]
        lon, lat = point_snapped.x, point_snapped.y
        ts = ds.sel(lon=lon, lat=lat, method='nearest').to_dataframe()
        ts['point'] = f"({lon}, {lat})"
        df = df.append(ts)

df.reset_index(inplace=True)