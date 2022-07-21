from shapely import geometry, ops
from shapely.geometry import Point, mapping, LineString, MultiPoint
import numpy as np
import rasterio
from rasterio import plot, mask
import networkx as nx
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from cartopy import crs
from scipy.ndimage import filters


def create_graph_and_paths(network_links):
    graph = nx.DiGraph()
    for index,row in network_links.iterrows():
        graph.add_edge(row['startnode'], row['endnode'], fid=row['fid'],length=row.length,dpn=row.dpn,
                       angle=row.angle,surface_cost=row.surface_cost,total_time=row.total_time)

    weighted_path_dpn = nx.dijkstra_path(graph, source=points[0], target=points[1], weight='dpn')
    weighted_path_surface_cost = nx.dijkstra_path(graph, source=points[0], target=points[1], weight='surface_cost')

    return weighted_path_dpn,weighted_path_surface_cost,graph


def create_path_gpd(weighted_path,graph):
    geom = []
    links = []
    length = []
    dpn = []
    angle = []
    surface_cost = []
    total_time = []
    first_node = weighted_path[0]
    for node in weighted_path[1:]:
        link_fid = graph.edges[first_node, node]['fid']
        links.append(link_fid)
        row = network_links.loc[network_links['fid'] == link_fid]
        geom.append(row['geometry'].cascaded_union)
        length.append(row.length.values[0])
        dpn.append(row.dpn.values[0])
        angle.append(row.angle.values[0])
        surface_cost.append(row.surface_cost.values[0])
        total_time.append(row.total_time.values[0])
        first_node = node

    weighted_path_gpd = gpd.GeoDataFrame({'fid': links,'length':length,'dpn':dpn,'angle':angle,
                                          'surface_cost':surface_cost,'total_time':total_time, 'geometry': geom})
    return weighted_path_gpd

def smooth_linestring(path_gpd, smooth_sigma):
    geom = path_gpd['geometry'].tolist()
    multi_line = geometry.MultiLineString(geom)
    multi_line
    linestring = ops.linemerge(multi_line)
    linestring
    start_coord = linestring.coords[0]
    end_coord = linestring.coords[-1]
    smooth_x = np.array(filters.gaussian_filter1d(
        linestring.xy[0],
        smooth_sigma)
    )
    smooth_y = np.array(filters.gaussian_filter1d(
        linestring.xy[1],
        smooth_sigma)
    )
    smoothed_coords = np.hstack((smooth_x, smooth_y))
    smoothed_coords = zip(smooth_x, smooth_y)
    smoothed_coords= [i for i in smoothed_coords]
    smoothed_coords[0] = start_coord
    smoothed_coords[-1] = end_coord
    linestring_smoothed = LineString(smoothed_coords)
    new_path_gpd = gpd.GeoSeries({'geometry': linestring_smoothed})
    return new_path_gpd

def plot_network(background_map, study_area_shapely, network_nodes, network_links, weighted_path_forward,
                 weighted_path_backward):
    back_array = background_map.read(1)
    palette = np.array([value for key, value in background_map.colormap(1).items()])
    background_image = palette[back_array]
    bounds = background_map.bounds
    extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)

    fig = plt.figure(figsize=(3, 3), dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=crs.OSGB())

    # display background map
    ax.imshow(background_image, origin='upper', extent=extent, zorder=0)

    # displaying nodes
    network_nodes.plot(ax=ax, zorder=3, markersize=0.2, alpha=0.5)

    # displaying links
    network_links.plot(ax=ax, zorder=2, edgecolor='blue', linewidth=0.2, alpha=0.5)

    # display path
    weighted_path_forward.plot(ax=ax, zorder=4, edgecolor='red', linewidth=0.7, label='surface cost')
    weighted_path_backward.plot(ax=ax, zorder=5, edgecolor='orange', linewidth=0.7, label='surface cost')

    # set the extent to the study area
    # study_area_gpd.plot(ax=ax,zorder = 2)
    display_extent = ((study_area_shapely.bounds[0] - 100, study_area_shapely.bounds[2] + 100,
                       study_area_shapely.bounds[1] - 100, study_area_shapely.bounds[3] + 100))
    ax.set_extent(display_extent, crs=crs.OSGB())
    plt.show()


OS_National_Grids = gpd.read_file(
    os.path.join('OS-British-National-Grids-main', 'OS-British-National-Grids-main', 'os_bng_grids.gpkg'),
    layer='1km_grid')
study_area_shapely = OS_National_Grids[OS_National_Grids['tile_name'] == "SX7478"].geometry.cascaded_union

SX77_map = rasterio.open(
    os.path.join('OS Explorer Maps', 'Download_SX77-Haytor_2033809', 'raster-25k_4596071', 'sx', 'sx77.tif'))

# network_links = gpd.read_file(os.path.join('Study_area', 'SX7677','Final Networks', 'network_links_dpn.geojson'))
# network_nodes = gpd.read_file(os.path.join('Study_area','SX7677','Final Networks', 'network_nodes_dpn.geojson'))

network_links = gpd.read_file(os.path.join('Study_area', 'SX7478', 'network_links_dpn.geojson'))
network_nodes = gpd.read_file(os.path.join('Study_area','SX7478', 'network_nodes_dpn.geojson'))

# Points for Haytor

# from Haytor rocks to B3387
# points = ['int_dpn1', 'dpn_1127']

# from bottom_right to Haytor Down
# points = ['dpn_408', 'dpn_298']

# from top_righ to haytor rocks path
# points = ['al_1681','int_dpn1']

# Points for Hound Tor

# from Hound Tor to Holwall lawn
# points = ["int_dpn17", "dpn_271"]

# from Holwell Lawn to medieval settlement
points = ["al_659", "al_1017"]

weighted_path_dpn,weighted_path_surface_cost,graph = create_graph_and_paths(network_links)


weighted_path_dpn_gpd = create_path_gpd(weighted_path_dpn,graph)
weighted_path_surface_cost_gpd  = create_path_gpd(weighted_path_surface_cost ,graph)


linestring_smoothed_dpn = smooth_linestring(weighted_path_dpn_gpd, 2)
linestring_smoothed_surface_cost = smooth_linestring(weighted_path_surface_cost_gpd, 2)


plot_network(SX77_map, study_area_shapely, network_nodes, network_links, weighted_path_dpn_gpd,
             linestring_smoothed_dpn)
plot_network(SX77_map, study_area_shapely, network_nodes, network_links, weighted_path_surface_cost_gpd,
             linestring_smoothed_surface_cost)

# linestring_smoothed_dpn.to_file("Study_area/SX7478/linestring_smoothed_dpn'int_dpn17','dpn_271'.geojson", driver='GeoJSON',crs='EPSG:27700')
# linestring_smoothed_surface_cost.to_file("Study_area/SX7478/linestring_smoothed_surface_cost'int_dpn17','dpn_271'.geojson", driver='GeoJSON',crs='EPSG:27700')

linestring_smoothed_dpn.to_file("Study_area/SX7478/linestring_smoothed_dpn'al_659','al_1017'.geojson", driver='GeoJSON',crs='EPSG:27700')
linestring_smoothed_surface_cost.to_file("Study_area/SX7478/linestring_smoothed_surface_cost'al_659','al_1017'.geojson", driver='GeoJSON',crs='EPSG:27700')