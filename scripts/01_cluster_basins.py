
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import matplotlib.pyplot as plt
import logging
import os 
from pathlib import Path
import seaborn as sns
from matplotlib.colors import ListedColormap
import traceback

#== PARAMS ==#
try:
    # INPUT
    basins = snakemake.input.basin_geojson
    
    #PARAMS
    crs = snakemake.params.crs
    cluster_method = snakemake.params.method
    intersect_method = snakemake.params.touches
    plot = snakemake.params.plot
    savefig = snakemake.params.savefig
    test_list = snakemake.params.test_list
    fill_rings = snakemake.params.fill_rings
    
    #LOGGING
    logging.info('Snakemake object found, using params from snakemake...')
    logging.info('Plotting: %s', plot)
    logging.info('Saving figures: %s', savefig)
    logging.info('Test list: %s', test_list)
    logging.info('CRS: %s', crs)
    logging.info('Cluster method: %s', cluster_method)
    logging.info('Intersect method: %s', intersect_method)
    logging.info('Basins file: %s', basins)
    logging.info('Test: %s', test_list)
    logging.info('Fill rings: %s', fill_rings)
    
    
except:
    os.chdir('p:/moonshot2-casestudy/Wflow/africa')
    logging.info('No snakemake object found, using default params...')
    basins = r"data\2-interim\GIS\basins_mainland_and_madagascar.geojson"
    cluster_method = "domain_method" #{domain_method, ...}
    intersect_method = "centroid" #{centroid, majority, ...}
    savefig = True
    crs = "EPSG:4326"
    plot = False
    test = False
    fill_rings = False
    test_list = None
    
#create logger
logger = logging.getLogger(__name__)
logging.basicConfig(filename=Path(os.path.join('data', '0-log', f'cluster_basins.log')).as_posix(), filemode='w', level=logging.INFO)
logging.info('Working directory: %s', os.getcwd())

interim_dir = Path(r'data\2-interim').as_posix()

#useful data 
africa = gpd.read_file(Path(r"data\1-external\administrative\Africa.shp").as_posix(), driver='ESRI Shapefile')

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
            return basins.sort_values(by=col, ascending=False)
        
        except:
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

# cluster basins
def cluster_basins(sorted_basin_data, mapped_bbox, filled_data, intersect_method, plot=False):
            filled_data['cluster_id'] = None

            clustered = set()
            cluster_dict = {}

            for n, i in enumerate(sorted_basin_data.index):
                try:
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


                    # Update 'cluster_id' in both 'sorted_basin_data' and 'filled_data'
                    sorted_basin_data.loc[cluster, 'cluster_id'] = i
                    filled_data.loc[cluster, 'cluster_id'] = i

                    if plot==True:
                        plot_basins(basin, bbox, filled_data.loc[cluster])

                    clustered.add(i)
                    clustered.update(cluster)
                    cluster_dict[i] = cluster
                except AttributeError as e:
                    logging.error('Error: %s', e)
                    continue

            sorted_basin_data['cluster_key'] = sorted_basin_data.index.to_series().apply(find_key, args=(cluster_dict,))
            
            return sorted_basin_data, cluster_dict, filled_data, cluster_dict

#function for test plotting
def plot_basins(basin: gpd.GeoDataFrame, #The basin that made the cluster
                bbox: gpd.GeoDataFrame, # 
                clustered: gpd.GeoDataFrame, 
                savefig:bool =False):
    
    plot_dir = r'data/5-visualization/cluster_plots'
    fig, ax = plt.subplots()
    gpd.GeoSeries(basin['geometry']).plot(ax=ax, color='blue', edgecolor='black')
    gpd.GeoSeries(bbox['geometry']).plot(ax=ax, edgecolor='red', facecolor='none')
    gpd.GeoSeries(clustered['geometry']).plot(ax=ax, color='green', edgecolor='black', alpha=0.5)
    gpd.GeoDataFrame(clustered['geometry']).dissolve().plot(ax=ax, color='none', edgecolor='pink', linewidth = 2, alpha=1)
    plt.title(f'Basin {basin.name} and intersecting basins')
    if savefig:
        os.makedirs(plot_dir, exist_ok=True)
        plt.savefig(os.path.join(plot_dir, f'cluster_{basin.name}.png'), dpi=400)
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

# lambda over the gdf to fill the column from the dict
def find_key(index, cluster_dict):
    for key, values in cluster_dict.items():
        if index in values:
            return key
    return None

def dissolve(gdf, by='cluster_id', test_list=None):
    if test_list is not None:
        gdf = gdf[gdf.cluster_id.isin(test_list)]
    return gdf.dissolve(by=by)

def to_polygon(geom):
    if geom.geom_type == 'Polygon':
        return geom
    elif geom.geom_type == 'MultiPolygon':
        # Convert the MultiPolygon to a list of Polygons
        polygons = [polygon for polygon in geom.geoms]
        # Return the Polygon with the maximum area
        return max(polygons, key=lambda part: part.area)
    else:
        # If it's a different type of geometry, raise an error
        raise ValueError(f'Cannot convert geometry of type {geom.geom_type} to Polygon')


# =============================================================================

#== MAIN ==#
if __name__ == "__main__":
    try:
        
        #read filedag
        basin_data = gpd.read_file(basins, crs=crs)
        logging.info('Basins data loaded: %s', basins)

        sorted_basin_data = sort_by_area(basin_data).set_index('fid')

        #map bbox to geometry column
        mapped_bbox = sorted_basin_data.bounds
        mapped_bbox['geometry'] = mapped_bbox.apply(lambda x: Polygon([(x.minx, x.miny), 
                                                                    (x.minx, x.maxy), 
                                                                    (x.maxx, x.maxy), 
                                                                    (x.maxx, x.miny)]), axis=1)
        #now a geodataframe of bboxes
        bbox_gdf = gpd.GeoDataFrame(mapped_bbox, crs=crs)

        #fill rings
        filled_data = sorted_basin_data.copy()

        #test fill rings
        if fill_rings == True:
            filled_data = fill_rings(sorted_basin_data)

        # Call the cluster_basins function
        sorted_basin_data, cluster_dict, filled_data, cluster_dict = cluster_basins(sorted_basin_data, 
                                                                                    mapped_bbox, 
                                                                                    filled_data, 
                                                                                    intersect_method, 
                                                                                    plot)

        if plot==True:
            # Create a colormap
            colors = sns.color_palette('Set3', n_colors=len(sorted_basin_data['cluster_key'].unique()))
            cmap = ListedColormap(colors)
            n_clusters = len(sorted_basin_data['cluster_key'].unique())
            fig, ax = plt.subplots()
            sorted_basin_data.plot(column='cluster_key', ax=ax, legend=True, cmap=cmap)
            ax.set_title(f'clustered by {cluster_method}, n = {n_clusters} clusters')

            if savefig:
                plt.savefig(Path(os.path.join(interim_dir, 'clustered_basins.png')).as_posix(), dpi=400)

        sorted_basin_data.to_file(Path(os.path.join(interim_dir, 'clustered_basins.geojson')).as_posix(), driver='GeoJSON')

        diss = dissolve(sorted_basin_data, by='cluster_key', test_list=test_list)
        diss = diss['geometry'].apply(to_polygon)

        diss.to_file(Path(os.path.join(interim_dir, 'dissolved_basins.geojson')).as_posix(), driver='GeoJSON')

        ls = list(diss.index.unique())
        
        with open(Path(os.path.join(interim_dir, 'cluster_list.txt')).as_posix(), 'w') as f:
            for item in ls:
                f.write("%s\n" % item)

        logging.info('Clustered basins saved to: %s', Path(os.path.join(interim_dir, 'clustered_basins.geojson')).as_posix())
        logging.info('Dissolved basins saved to: %s', Path(os.path.join(interim_dir, 'dissolved_basins.geojson')).as_posix())
    
    except Exception as e:
        logging.error('Error: %s', e)
        traceback.print_exc()
        exit(1)
    
