import traceback
import os
from pathlib import Path
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.geometry import box
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import sys
sys.path.append('c:/git/MOONSHOT_2_FOR_AFRICA')
from scripts.helper import setup_logging

# == FUNCTIONS ==#
# sort basins by area
def sort_by_area(basins: gpd.GeoDataFrame):
    try:
        return basins.sort_values(by="area_km2", ascending=False)
    except:
        try:
            cols = basins.columns
            col = [c for c in cols if "area" in c.lower()][0]
            return basins.sort_values(by=col, ascending=False)

        except:
            exit(1)

# function that returns rules based intersecting basins
def intersecting(bbox, targets, method="centroid"):
    '''
    For each bbox, if method:
    - centroid: find which other centroids fall within the bbox
    - ?? 
    '''
    if method == "centroid":
        # Convert to a projected CRS before calculating centroid
        targets = targets.to_crs("EPSG:3395").centroid
        targets = targets.to_crs(crs)
        intersecting_indices = targets.index[targets.intersects(bbox)]
    if method == "touches":
        intersecting_indices = targets.index[targets.touches(bbox)]
    
    return intersecting_indices

# cluster basins
def cluster_basins(logger, sorted_basin_data, intersect_method, recurse, plot, eliminate_solo):
    # map bbox to geometry column
    mapped_bbox = sorted_basin_data.bounds
    mapped_bbox["geometry"] = mapped_bbox.apply(
        lambda x: Polygon(
            [(x.minx, x.miny), (x.minx, x.maxy), (x.maxx, x.maxy), (x.maxx, x.miny)]
        ),
        axis=1,
    )
    clustered_basins = gpd.GeoDataFrame(
        index=sorted_basin_data.index,
        columns=["cluster_key"],
        crs=sorted_basin_data.crs,
        geometry=sorted_basin_data.geometry,
    )
    
    clustered = set()
    cluster_dict = {}

    for n, i in enumerate(sorted_basin_data.index):
        if i in clustered:
            continue

        basin = sorted_basin_data.loc[i]
        bbox = mapped_bbox.loc[i]

        #find if any neighboring centroids are in the bbox
        cluster = intersecting(bbox["geometry"], sorted_basin_data, intersect_method)

        sorted_basin_data.loc[cluster, "cluster_id"] = i
 
        
        clustered.add(i)
        clustered.update(cluster)
        cluster_dict[i] = cluster
        # logger.info("Cluster %s: %s", i, cluster)


    clustered_basins["cluster_key"] = sorted_basin_data.index.to_series().apply(
        find_key, args=(cluster_dict,)
    )

    # Dissolve to simple polygons. Sort the values so the output text file is ordered by area
    dissolved_basins = dissolve(clustered_basins, by="cluster_key", test_list=test_list)
    clustered_basins["cluster_length"] = clustered_basins["cluster_key"].apply(
        lambda x: len(cluster_dict[x])
    )
    #clusters of len 1
    len_1 = clustered_basins[clustered_basins["cluster_length"] == 1]
    logger.info("Clusters of length 1: %s", len(len_1.index))
    
    if recurse is not None and recurse > 0:
        dissolved_basins, clustered_basins, cluster_dict = recursive_cluster(
            dissolved_basins, clustered_basins, cluster_dict, logger, recurse
        )
    else:
        logger.info("No recursion is enforced, clustering only what falls within the bounds of the largest basins")

    return clustered_basins, dissolved_basins, cluster_dict

def plot_diss(diss, clustered, mapped_bbox, savefig=False):
    fig, ax = plt.subplots()
    for i in diss.index:
        basin = diss.loc[i]
        bbox = mapped_bbox.loc[i]
        gpd.GeoSeries(basin["geometry"]).plot(ax=ax, color="blue", edgecolor="black")
        gpd.GeoSeries(bbox["geometry"]).plot(ax=ax, edgecolor="red", facecolor="none")
        gpd.GeoSeries(clustered["geometry"]).plot(
            ax=ax, color="green", edgecolor="black", alpha=0.5
        )
        gpd.GeoDataFrame(clustered["geometry"]).dissolve().plot(
            ax=ax, color="none", edgecolor="pink", linewidth=2, alpha=1
        )
        plt.title(f"Basin {i} and intersecting basins")
        if savefig:
            plt.savefig(os.path.join(plot_dir, f"cluster_{i}.png"), dpi=400)
        plt.show()

def recursive_cluster(diss, clustered, cluster_dict, logger, recurse, intersect_method="centroid", plot=True):
    '''Accomplish More Clustering by:
            1. Using the large clusters buffered to see if they intersect more basins eliminating small clusters
    '''
    
    logger.info("Recursively clustering basins")

    new_diss = None
    for _r in range(recurse + 1):
        logger.info("Recursion %s", _r)
        if new_diss is None:
            new_diss = diss
        else:
            new_diss = new_diss
        
        mapped_bbox = new_diss.bounds
        
        mapped_bbox["geometry"] = mapped_bbox.apply(
            lambda x: Polygon(
                [(x.minx, x.miny), (x.minx, x.maxy), (x.maxx, x.maxy), (x.maxx, x.miny)]
            ),
            axis=1,
        )
        
        for i in diss.index:
            logger.info("Cluster %s", i)
            basin = diss.loc[i]
            bbox = mapped_bbox.loc[i]
            cluster = intersecting(bbox["geometry"], new_diss, intersect_method)
            logger.info("Cluster %s: %s", i, cluster)
            new_adds = []
            for c in cluster:
                new_adds.extend(cluster_dict[c])

            diss.loc[cluster, "cluster_id"] = i
            cluster_dict[i] = list(cluster_dict[i]) + list(new_adds)
            
            for c in cluster:
                for id in cluster_dict[c]:
                    clustered.loc[id, "cluster_key"] = i
                cluster_dict.pop(c)
        if plot:
            plot_diss(diss, clustered, mapped_bbox)
        
        
        
# function for test plotting
def plot_basins(
    basin: gpd.GeoDataFrame,  # The basin that made the cluster
    bbox: gpd.GeoDataFrame,  #
    clustered: gpd.GeoDataFrame,
    savefig: bool = False,
):

    plot_dir = r"data/5-visualization/cluster_plots"
    fig, ax = plt.subplots()
    gpd.GeoSeries(basin["geometry"]).plot(ax=ax, color="blue", edgecolor="black")
    gpd.GeoSeries(bbox["geometry"]).plot(ax=ax, edgecolor="red", facecolor="none")
    gpd.GeoSeries(clustered["geometry"]).plot(
        ax=ax, color="green", edgecolor="black", alpha=0.5
    )
    gpd.GeoDataFrame(clustered["geometry"]).dissolve().plot(
        ax=ax, color="none", edgecolor="pink", linewidth=2, alpha=1
    )
    plt.title(f"Basin {basin.name} and intersecting basins")
    if savefig:
        os.makedirs(plot_dir, exist_ok=True)
        plt.savefig(os.path.join(plot_dir, f"cluster_{basin.name}.png"), dpi=400)
    plt.show()


def area_threshold_merge(minimum_area, gdf):
    # Merge clusters below threshold into their neighbours
    threshold_area = minimum_area * float(gdf["area_km2"].max())

    # Separate the small clusters from diss based on threshold area
    small_polygons = gdf[gdf["area_km2"] <= threshold_area]
    gdf = gdf[gdf["area_km2"] > threshold_area]

    # Loop through small clusters, find id of closest neighbour and merge
    for _, small_polygon in small_polygons.iterrows():
        distances = gdf.geometry.apply(
            lambda geom: small_polygon.geometry.distance(geom)
        )
        id_closest = distances.idxmin()
        gdf.at[id_closest, "geometry"] = gdf.at[id_closest, "geometry"].union(
            small_polygon["geometry"]
        )
    return gdf

# lambda over the gdf to fill the column from the dict
def find_key(index, cluster_dict):
    for key, values in cluster_dict.items():
        if index in values:
            return key
    return None

def dissolve(gdf, by="cluster_id", test_list=None):
    if test_list is not None:
        gdf = gdf[gdf.cluster_id.isin(test_list)]
    gdf = gdf.dissolve(by=by)
    gdf["geometry"] = gdf["geometry"].apply(to_polygon)
    gdf = gdf.to_crs(crs)
    gdf["area_km2"] = gdf.area
    gdf = gdf.sort_values(by="area_km2", ascending=False)
    return gdf

def to_polygon(geom):
    if geom.geom_type == "Polygon":
        return geom
    elif geom.geom_type == "MultiPolygon":
        # Convert the MultiPolygon to a list of Polygons
        polygons = [polygon for polygon in geom.geoms]
        # Return the Polygon with the maximum area
        return max(polygons, key=lambda part: part.area)
    else:
        # If it's a different type of geometry, raise an error
        raise ValueError(f"Cannot convert geometry of type {geom.geom_type} to Polygon")

def main(
    logger,
    basins: str|Path, 
    crs: str = "EPSG:4326", 
    cluster_method:str = "domain_method", 
    intersect_method: str = "centroid", 
    plot: bool = False, 
    savefig: bool = False, 
    test_list: list = None, 
    recursive: int = 0, 
    minimum_area: float = None,
    fill_rings: bool = False,
    eliminate_solo: bool = False):
    
    basin_data = gpd.read_file(basins, crs=crs)
    logger.info("Basins data loaded: %s", basins)

    sorted_basin_data = sort_by_area(basin_data).set_index("fid")

    # Call the cluster_basins function
    sorted_basin_data, diss, cluster_dict= cluster_basins(
        logger, sorted_basin_data, intersect_method, recursive, plot, eliminate_solo
    )

    if test_list:
        sorted_basin_data = sorted_basin_data[
            sorted_basin_data["cluster_id"].isin(test_list)
        ]

    if plot:
        # Create a colormap
        colors = sns.color_palette(
            "Set3", n_colors=len(sorted_basin_data["cluster_key"].unique())
        )
        cmap = ListedColormap(colors)
        n_clusters = len(sorted_basin_data["cluster_key"].unique())
        fig, ax = plt.subplots()
        sorted_basin_data.plot(column="cluster_key", ax=ax, legend=True, cmap=cmap)
        ax.set_title(f"clustered by {cluster_method}, n = {n_clusters} clusters")

        if savefig:
            plt.savefig(
                Path(
                    os.path.join(visualization_dir, "clustered_basins.png")
                ).as_posix(),
                dpi=400,
            )

    sorted_basin_data.to_file(
        Path(os.path.join("data/2-interim", "clustered_basins.geojson")).as_posix(),
        driver="GeoJSON",
    )

    if minimum_area:
        diss = area_threshold_merge(minimum_area, diss)

    diss.to_file(
        Path(os.path.join("data/2-interim", "dissolved_basins.geojson")).as_posix(),
        driver="GeoJSON",
    )

    ls = list(diss.index.unique())

    with open(
        Path(os.path.join("data/2-interim", "cluster_list.txt")).as_posix(), "w"
    ) as f:
        for item in ls:
            if item not in [1893, 1892, 1025, 1784, 1031, 1027]:
                f.write("%s\n" % item)
                os.makedirs(
                    os.path.join(output_dir, "models", f"item"), exist_ok=True
                )

    logger.info(
        "Clustered basins saved to: %s",
        Path(os.path.join("data/2-interim", "clustered_basins.geojson")).as_posix(),
    )
    logger.info(
        "Dissolved basins saved to: %s",
        Path(os.path.join("data/2-interim", "dissolved_basins.geojson")).as_posix(),
    )

# == MAIN ==#
if __name__ == "__main__":
    logger = setup_logging("data/0-log", "cluster_basins.log")
    try:
        if 'snakemake' in globals():
            smk = globals()['snakemake']
            logger.info("Snakemake object found, using params from snakemake...")
            main(
                logger=logger, 
                basins=smk.input.basin_geojson, 
                crs=smk.params.crs, 
                method=smk.params.method, 
                touches=smk.params.touches, 
                plot=smk.params.plot, 
                savefig=smk.params.savefig, 
                test_list=smk.params.test_list, 
                recursive=smk.params.recursive, 
                fill_rings=smk.params.fill_rings, 
                eliminate_solo=smk.params.eliminate_solo
            )
        else: 
            logger.info("No snakemake object found, using default params...")
            os.chdir('p:/moonshot2-casestudy/Wflow/africa')
            cluster_method = "domain_method"  # {domain_method, ...}
            intersect_method = "centroid"  # {centroid, majority, ...}
            savefig = False
            crs = "EPSG:4326"
            plot = True
            test = True
            fill_rings = False
            test_list = None
            minimum_area = None  # 0.001
            basins=Path("data", "2-interim", "GIS", "basins_mainland_and_madagascar.geojson")
            recursive = 3
            eliminate_solo = True
            main(
                logger=logger, 
                basins=basins, 
                crs=crs, 
                cluster_method=cluster_method, 
                intersect_method=intersect_method, 
                plot=plot, 
                savefig=savefig, 
                test_list=test_list, 
                fill_rings=fill_rings, 
                minimum_area=minimum_area, 
                recursive=recursive,
                eliminate_solo=eliminate_solo
            )

    except Exception as e:
        logger.error("Error: %s", e)
        logger.error(traceback.format_exc())
        exit(1)
