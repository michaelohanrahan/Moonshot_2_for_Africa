
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import matplotlib.pyplot as plt
import logging
import os 
import seaborn as sns
from matplotlib.colors import ListedColormap

#create logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

os.chdir(r'p:\moonshot2-casestudy\Wflow\africa')
logging.info('Working directory: %s', os.getcwd())

#== PARAMS ==#
try:
    basins = snakemake.input.basin_geojson
    crs = snakemake.params.crs
    cluster_method = snakemake.params.method
    plot = snakemake.params.plot
    logging.info('Snakemake object found, using params from snakemake...')
    logging.info('Plot: %s', plot)
    
except:
    logging.info('No snakemake object found, using default params...')
    basins = r"data\2-interim\GIS\basins_mainland_and_madagascar.geojson"
    cluster_method = "domain_method" #{domain_method, ...}
    intersect_method = "centroid" #{centroid, majority, ...}
    crs = "EPSG:4326"
    plot = True

interim_dir = r'data\2-interim'

#useful data 
africa = gpd.read_file(r"data\1-external\administrative\Africa.shp", driver='ESRI Shapefile')

#log input params
logging.info('Basins file: %s', basins)
logging.info('CRS: %s', crs)
logging.info('Cluster method: %s', cluster_method)
logging.info('Intersect method: %s', intersect_method)

# =============================================================================

#== FUNCTIONS ==#
#sort basins by area
def sort_by_area(basins: gpd.GeoDataFrame):
    try:
        
        return basins.sort_values(by='area_km2', ascending=False)
    except:
        try:
            cols = basins.columns
            col = [c for c in cols if 'area' in c.lower()][0]
            print('col "area_km2" not found, using col:', col, 'instead.')
            return basins.sort_values(by=col, ascending=False)
        
        except:
            print('Error: no area column found in basins dataframe. Exiting...')
            exit(1)

#function that retruns rules based intersecting basins
def intersecting(bbox, targets, method='centroid'):
    if method == 'centroid':
        # Convert to a projected CRS before calculating centroid
        targets = targets.to_crs('EPSG:3395').centroid

    targets = targets.to_crs(crs)
    intersecting_indices = targets.index[targets.intersects(bbox)]
    logging.info('Indices of basins intersecting bbox: %s', intersecting_indices)
    return intersecting_indices

#function for test plotting
def plot_basins(basin: gpd.GeoDataFrame, #The basin that made the cluster
                bbox: gpd.GeoDataFrame, # 
                clustered: list):
    fig, ax = plt.subplots()
    gpd.GeoSeries(basin['geometry']).plot(ax=ax, color='blue', edgecolor='black')
    gpd.GeoSeries(bbox['geometry']).plot(ax=ax, edgecolor='red', facecolor='none')
    gpd.GeoSeries(clustered['geometry']).plot(ax=ax, color='green', edgecolor='black', alpha=0.5)
    
    plt.show()

#function for filling rings in the basins dataset
def fill_rings(basins):
    '''
    Iterate the basins. If the geometry is a Polygon, extract the interior rings and create new polygons.
    If the geometry is a MultiPolygon, extract the interior rings and create new polygons.
    
    Use these new polygons to fill the interior rings in the original basins dataset.
    
    '''
    
    new_polygons = []
    
    for idx, row in basins.iterrows():
        if isinstance(row['geometry'], Polygon):
            # Extract rings and create new polygons
            for interior in row['geometry'].interiors:
                new_polygons.append(Polygon(interior))

            # Fill rings for the Polygon
            basins.loc[idx, 'geometry'] = row['geometry'].buffer(0)
        
        elif isinstance(row['geometry'], MultiPolygon):
            # Handle MultiPolygon
            updated_polygons = []
            for polygon in row['geometry'].geoms:
                # Extract rings and create new polygons
                for interior in polygon.interiors:
                    # if interior.area > 1:
                    new_polygons.append(Polygon(interior))

                # Fill rings for the Polygon and add to the list
                if polygon.area > 1:
                    updated_polygons.append(polygon.buffer(0))

            # Update the MultiPolygon with the filled polygons
            basins.loc[idx, 'geometry'] = MultiPolygon(updated_polygons)
            
    print('new polygons:', new_polygons)
    fig, ax = plt.subplots(figsize=(10, 10))
    gpd.GeoDataFrame(africa, columns=['geometry']).plot(ax=ax, edgecolor='black', facecolor='grey', alpha=0.7)
    gpd.GeoDataFrame(new_polygons, columns=['geometry']).plot(ax=ax, color='red', edgecolor='white', alpha=0.5)
    plt.title('Endorheic Basins of Mainland African Continent and Madagascar', fontsize=14)
    plt.suptitle(f'n = {len(new_polygons)} polygons created from filling rings in the basins dataset.', y=0.88, fontsize=12)
    plt.savefig(r'data/5-visualization/filled_rings.png', dpi=400)
    plt.show()
    
    # Remove entries that are within the new polygons
    for new_polygon in new_polygons:
        for idx, row in basins.iterrows():
            if row['geometry'].within(new_polygon):
                basins = basins.drop(idx)

    return basins

# =============================================================================

#== MAIN ==#
# if __name__ == "__main__":
basin_data = gpd.read_file(basins, crs=crs)
sorted_basin_data = sort_by_area(basin_data).set_index('fid')

#map bbox to geometry column
mapped_bbox = sorted_basin_data.bounds
mapped_bbox['geometry'] = mapped_bbox.apply(lambda x: Polygon([(x.minx, x.miny), 
                                                               (x.minx, x.maxy), 
                                                               (x.maxx, x.maxy), 
                                                               (x.maxx, x.miny)]), axis=1)
#now a geodataframe of bboxes
bbox_gdf = gpd.GeoDataFrame(mapped_bbox, crs=crs)

#test fill rings
filled_data = sorted_basin_data.copy()
# filled_data = fill_rings(sorted_basin_data)

filled_data['cluster_id'] = None

clustered = set()
cluster_dict = {}
count = 0

for n, i in enumerate(sorted_basin_data.index):
    if i in clustered:
        continue

    # Mask the gdf by the growing cluster set
    if len(clustered) > 0:
        filled_data = filled_data.loc[~filled_data.index.isin(clustered)] 

    logging.info('Basin: %s', i)
    basin = sorted_basin_data.loc[i]
    bbox = mapped_bbox.loc[i]

    cluster = intersecting(bbox['geometry'], filled_data, intersect_method)

    logging.info('(i, len(cluster)) %s', (i,len(cluster)))

    plot_basins(basin, bbox, filled_data.loc[cluster])

    # Update 'cluster_id' in both 'sorted_basin_data' and 'filled_data'
    sorted_basin_data.loc[cluster, 'cluster_id'] = i
    filled_data.loc[cluster, 'cluster_id'] = i
    
    clustered.add(i)
    clustered.update(cluster)
    cluster_dict[i] = cluster
    print(clustered)
    count += 1

# lambda over the gdf to fil lthe columnf from the dict
def find_key(index):
    for key, values in cluster_dict.items():
        if index in values:
            return key
    return None

sorted_basin_data['cluster_key'] = sorted_basin_data.index.to_series().apply(find_key)

if plot:
    # Create a colormap
    colors = sns.color_palette('Set3', n_colors=len(sorted_basin_data['cluster_key'].unique()))
    cmap = ListedColormap(colors)

    fig, ax = plt.subplots()
    sorted_basin_data.plot(column='cluster_key', ax=ax, legend=True, cmap=cmap)
    plt.show()

sorted_basin_data.to_file(os.path.join(interim_dir, 'clustered_basins.geojson'), driver='GeoJSON')
logging.info('Clustered basins saved to: %s', os.path.join(interim_dir, 'clustered_basins.geojson'))