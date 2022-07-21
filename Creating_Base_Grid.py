
# packages needed
import shapely.affinity
from shapely.geometry import Point,mapping, LineString, Polygon
import numpy as np
import rasterio
from rasterio import plot, mask
from rtree import index
import networkx as nx
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from cartopy import crs
from math import atan, degrees

# functions

def create_nodes(study_area_shapely,elevation_mask,transform_index):
    count_x = np.arange(study_area_shapely.bounds[0], study_area_shapely.bounds[2] + 25, 25, dtype=int)
    count_y = np.arange(study_area_shapely.bounds[1], study_area_shapely.bounds[3] + 25, 25, dtype=int)

    points = []
    height = []
    for i in count_x:
        for j in count_y:
            points.append(Point(i,j))
            point_index = rasterio.transform.rowcol(transform_index, i,
                                                    j)
            point_elevation = elevation_mask[0][point_index]
            height.append(point_elevation)
    nodes_gpd = gpd.GeoSeries(points)

    G = nodes_gpd.geometry.apply(lambda geom: geom.wkb)
    nodes_gpd =nodes_gpd.loc[G.drop_duplicates().index]

    network_nodes = gpd.GeoDataFrame(geometry=nodes_gpd)
    network_nodes['fid'] = range(1,len(network_nodes)+1)
    network_nodes['fid'] = 'al_' + network_nodes['fid'].astype(str)
    network_nodes['height'] = height
    return network_nodes


def create_links(network_nodes):
    lines = []
    start_node = []
    end_node = []
    length = []
    angles = []
    climb_time = []
    for index,row in network_nodes.iterrows():
        buffer = row.geometry.buffer(37)
        intersection = network_nodes.intersection(buffer)
        intersection = intersection[~intersection.is_empty]
        for j in intersection[1:]:
            # get start and end node
            start_node.append(row.fid)
            end_point = network_nodes.loc[network_nodes['geometry']== j,'fid'].iloc[0]
            end_node.append(end_point)
            start_node.append(end_point)
            end_node.append(row.fid)

            # get the line geometry
            line_forwards = LineString([row.geometry,j])
            #line_backwards = LineString([j,row.geometry])
            lines.append(line_forwards)
            lines.append(line_forwards)

            # get the length of the line
            length.append(line_forwards.length)
            length.append(line_forwards.length)

            # get the angle of elevation change and climb time
            start_height_forward = end_height_backward = row.height
            end_height_forward = start_height_backward = network_nodes.loc[network_nodes['geometry']== j,'height'].iloc[0]

            if start_height_forward > end_height_forward:
                change_height = start_height_forward-end_height_forward
                climb_time.append(0)
            else:
                change_height = end_height_forward-start_height_forward
                climb_time.append(change_height/10)

            angle = degrees(atan((change_height/line_forwards.length)))
            angles.append(angle)

            if start_height_backward > end_height_backward:
                change_height = start_height_backward-end_height_backward
                climb_time.append(0)
            else:
                change_height = end_height_backward-start_height_backward
                climb_time.append(change_height/10)

            angle = degrees(atan((change_height/line_forwards.length)))
            angles.append(angle)

    links_fid = range(1,len(lines)+1)
    network_links = gpd.GeoDataFrame({'fid':links_fid,'startnode':start_node,
                                      'endnode':end_node,'length':length,'angle':angles,
                                      'climb_time_forward':climb_time,'geometry':lines})

    G = network_links['geometry'].apply(lambda geom: geom.wkb)
    network_links=network_links.loc[G.drop_duplicates().index]
    network_links = network_links.drop(network_links[(network_links.geometry.length == 0 )].index)

    return network_links

def plot_network(background_map,study_area_shapely,elevation_mask,transform_index,network_nodes,network_links):
    back_array = background_map.read(1)
    palette = np.array([value for key, value in background_map.colormap(1).items()])
    background_image = palette[back_array]
    bounds = background_map.bounds
    extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
    fig = plt.figure(figsize=(3, 3), dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=crs.OSGB())

    # display background map
    ax.imshow(background_image, origin='upper', extent=extent, zorder=0)

    # display elevation
    rasterio.plot.show(elevation_mask, alpha=0.6, transform=transform_index, ax=ax, zorder=1,
                       cmap='terrain')

    # displaying nodes
    network_nodes.plot(ax=ax,zorder = 3,markersize=0.2)

    #displaying links
    network_links.plot(ax=ax,zorder = 2,edgecolor='blue', linewidth=0.2)

    #set the extent to the study area
    display_extent = ((study_area_shapely.bounds[0]- 100, study_area_shapely.bounds[2]+ 100,
                       study_area_shapely.bounds[1]- 100, study_area_shapely.bounds[3]+ 100))

    ax.set_extent(display_extent, crs=crs.OSGB())
    plt.show()

#files to import

OS_National_Grids = gpd.read_file(os.path.join('OS-British-National-Grids-main','OS-British-National-Grids-main','os_bng_grids.gpkg'),layer='1km_grid')

study_area_shapely = OS_National_Grids[OS_National_Grids['tile_name'] == "SX7478"].geometry.cascaded_union

SX77_map = rasterio.open(
    os.path.join('OS Explorer Maps', 'Download_SX77-Haytor_2033809', 'raster-25k_4596071', 'sx', 'sx77.tif'))

elevation = rasterio.open(os.path.join('OS Elevation','SX77_elevation','terrain-5-dtm_4616587','sx','sx77ne_nw','w001001.adf'))

def limit_elevation(elevation,study_area_shapely):
    elevation_study_area = shapely.affinity.scale(study_area_shapely,xfact=1.1, yfact=1.1,origin='center')
    study_area = mapping(elevation_study_area)
    elevation_mask, transform_index = mask.mask(elevation,[study_area], filled=False, crop=False)
    return elevation_mask,transform_index

# # creating elevation mask of the study area
# elevation_study_area = shapely.affinity.scale(study_area_shapely,xfact=1.1, yfact=1.1,origin='center')
# study_area = mapping(elevation_study_area)
# elevation_mask, transform_index = mask.mask(elevation,[study_area], filled=False, crop=False)

elevation_mask,transform_index = limit_elevation(elevation,study_area_shapely)

network_nodes = create_nodes(study_area_shapely,elevation_mask,transform_index)

network_links = create_links(network_nodes)

plot_network(SX77_map,study_area_shapely,elevation_mask,transform_index,network_nodes,network_links)

# network_links.to_file("Study_area/SX7677/Final Networks/network_links_al.geojson", driver='GeoJSON',crs='EPSG:27700')
# network_nodes.to_file("Study_area/SX7677/Final Networks/network_nodes_al.geojson", driver='GeoJSON',crs='EPSG:27700')

network_links.to_file("Study_area/SX7478/network_links_al.geojson", driver='GeoJSON',crs='EPSG:27700')
network_nodes.to_file("Study_area/SX7478/network_nodes_al.geojson", driver='GeoJSON',crs='EPSG:27700')

