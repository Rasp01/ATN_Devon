from shapely.geometry import Point, LineString
import numpy as np
import rasterio
from rtree import index
import networkx as nx
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from cartopy import crs
import xml.etree.ElementTree as ET
import gpxpy.gpx
from pyproj import CRS
from pyproj import Transformer


def create_network(gml, gml_gpd):
    root = gml.getroot()
    startnodes = [startNode.attrib for startNode in
                  root.iter('{http://namespaces.os.uk/networks/detailedPathNetwork/1.0}startNode')]
    gml_gpd['startNodes']  = [nodes['{http://www.w3.org/1999/xlink}href'] for nodes in startnodes]
    endnodes = [endNode.attrib for endNode in
                root.iter('{http://namespaces.os.uk/networks/detailedPathNetwork/1.0}endNode')]
    gml_gpd['endNodes'] = [nodes['{http://www.w3.org/1999/xlink}href'] for nodes in endnodes]

    nodes_id = [nodes.attrib for nodes in
                root.iter('{http://namespaces.os.uk/networks/detailedPathNetwork/1.0}RouteNode')]
    node_ids= [node['{http://www.opengis.net/gml/3.2}id'] for node in nodes_id]

    network_coords = [geom.text for geom in root.iter('{http://www.opengis.net/gml/3.2}pos')]
    node_coordinates = []
    for i in network_coords:
        coords = i.split(" ")
        coords = [float(j) for j in coords]
        node_coordinates.append(coords)

    node_geom = [Point(node) for node in node_coordinates]

    network_nodes = gpd.GeoDataFrame({'fid': nodes_id, 'geometry': node_geom})

    graph = nx.Graph()
    for index, row in gml_gpd.iterrows():
        graph.add_edge(row['startNodes'], row['endNodes'], fid=row['gml_id'], weight=row['planimetricLength'])

    return graph, network_nodes, gml_gpd, node_ids, node_coordinates


def import_gpx(gpx_file):
    gpx = gpxpy.parse(gpx_file)

    wsg84 = CRS.from_epsg(4326)
    osgb36 = CRS.from_epsg(27700)
    transformer = Transformer.from_crs(wsg84, osgb36)

    points = []
    if gpx.waypoints:
        waypoints = gpx.waypoints
    else:
        routes = gpx.routes
        waypoints = routes[0].points

    for waypoint in waypoints:
        points.append(transformer.transform(waypoint.latitude, waypoint.longitude))

    route_gpd = gpd.GeoDataFrame(index=[0], crs='epsg:27700', geometry=[LineString(points)])
    #route_gpd.plot()

    bounds = []
    for bound in gpx.bounds:
        bounds.append(bound)
    top_right = (transformer.transform(bounds[1], bounds[3]))
    bottom_left = (transformer.transform(bounds[0], bounds[2]))

    return route_gpd, points, top_right, bottom_left


def get_nearest_nodes(node_ids, node_coordinates, first_point, point, global_nodes,global_missing_links):
    idx = index.Index()

    # set the bounds for the index
    for i in range(len(node_ids)):
        left, bottom, right, top = (node_coordinates[i][0], node_coordinates[i][1],
                                    node_coordinates[i][0], node_coordinates[i][1])
        idx.insert(i, (left, bottom, right, top))

    for i in idx.nearest(first_point, 1):
        first_coordinate = node_coordinates[i]
        # need to add # so that it is identified in with the link
        first_node = ("#" + node_ids[i])
        global_nodes.append(Point(first_coordinate))
    for i in idx.nearest(point, 1):
        last_coordinate = node_coordinates[i]
        # need to add # so that it is identified in with the link
        last_node = ("#" + node_ids[i])
        global_nodes.append(Point(last_coordinate))

    if first_coordinate == last_coordinate:
        print('Same node')
        #plt.plot(LineString(zip(first_point, point)))
        global_missing_links.append(LineString([first_point, point]))
    return first_coordinate, first_node, last_coordinate, last_node


def gpx_to_path(first_node, last_node, g, path_network, global_geom, global_links):
    path = nx.dijkstra_path(g, first_node, last_node, weight='weight')
    if len(path) == 1:
         print('path list only contains 1 node')
    else:
        geom = []
        links = []

        first_node = path[0]
        for node in path[1:]:
            link_fid = g.edges[first_node, node]['fid']
            links.append(link_fid)
            global_links.append(link_fid)
            row = path_network.loc[path_network['gml_id'] == link_fid]
            geom.append(row['geometry'].cascaded_union)
            global_geom.append(row['geometry'].cascaded_union)
            first_node = node

        path_gpd = gpd.GeoDataFrame({'fid': links, 'geometry': geom})
        #print('Plotting Path')
        #path_gpd.plot()
        return path_gpd


# def plot_map(sheepstor_map, first_coordinate, last_coordinate, coords,
#              bluebell_walk, path_network, path_nodes, path_gpd):
#     back_array = sheepstor_map.read(1)
#     palette = np.array([value for key, value in sheepstor_map.colormap(1).items()])
#     background_image = palette[back_array]
#     bounds = sheepstor_map.bounds
#     extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
#     fig = plt.figure(figsize=(3, 3), dpi=500)
#     ax = fig.add_subplot(1, 1, 1, projection=crs.OSGB())
#     ax.imshow(background_image, origin='upper', extent=extent, zorder=0)
#     plt.scatter(first_coordinate[0], first_coordinate[1], color='green', s=1, zorder=5)
#     plt.scatter(last_coordinate[0], last_coordinate[1], color='green', s=1, zorder=5)
#     bluebell_walk.plot(ax=ax, edgecolor='red', linewidth=0.5, zorder=2)
#     plt.scatter(*zip(*coords), color='red', s=1, zorder=3)
#     path_network.plot(ax=ax, edgecolor='blue', linewidth=0.5, zorder=4)
#     path_nodes.plot(ax=ax, color='blue', markersize=1, zorder=4)
#     path_gpd.plot(ax=ax, edgecolor='green', linewidth=0.5, zorder=5)
#     display_extent = ((first_coordinate[0] - 500, first_coordinate[0] + 500,
#                        first_coordinate[1] - 500, first_coordinate[1] + 500))
#     ax.set_extent(display_extent, crs=crs.OSGB())
#     plt.show()

def calculate_missing_sections(global_missing_links, path_gpd):
    missing_links = []
    path_geoseries = gpd.GeoSeries(path_gpd['geometry'])
    # account for missing nodes that are on the path due to being so small there nodes were the same
    for link in global_missing_links:
        if any(path_geoseries.intersects(link)):
            print("Line is present on final path ")
        else:
            missing_links.append(link)
        # if path_geoseries.intersects(link) == True:
        #     print("Line is present on final path ")
        # else:
        #     missing_links.append(link)
    #global_missing_links_gpd = gpd.GeoSeries(global_missing_links)
    global_missing_links_gpd = gpd.GeoSeries(missing_links)

    return global_missing_links_gpd

def plot_missing_links(global_missing_links_gpd, global_missing_links):
    global_missing_links_gpd.plot()

def plot_global_map(raster_map, global_nodes_gpd, points,
             route_gpd, path_network, path_nodes, path_gpd, top_right, bottom_left, global_missing_links_gpd):
    back_array = raster_map.read(1)
    palette = np.array([value for key, value in raster_map.colormap(1).items()])
    background_image = palette[back_array]
    bounds =  raster_map.bounds
    extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
    fig = plt.figure(figsize=(3, 3), dpi=500)

    ax = fig.add_subplot(1, 1, 1, projection=crs.OSGB())
    ax.imshow(background_image, origin='upper', extent=extent, zorder=0)
    route_gpd.plot(ax=ax, edgecolor='blue', linewidth=0.5, zorder=2)
    plt.scatter(*zip(*points), color='blue', s=1, zorder=3)

    # path_network.plot(ax=ax, edgecolor='blue', linewidth=0.5, zorder=4)
    # path_nodes.plot(ax=ax, color='blue', markersize=1, zorder=4)

    global_missing_links_gpd.plot(ax=ax, color='red', linewidth=1, zorder=6)

    path_gpd.plot(ax=ax, edgecolor='green', linewidth=1, zorder=5)
    global_nodes_gpd.plot(ax=ax, color='yellow', markersize=1, zorder=5)

    display_extent = ((bottom_left[0] - 100, top_right[0] + 100,
                       bottom_left[1] - 100, top_right[1] + 100))
    ax.set_extent(display_extent, crs=crs.OSGB())
    plt.show()


def main():
    sheepstor_map = rasterio.open(
        os.path.join('OS Explorer Maps', 'Download_South+Dartmoor_2004150', 'raster-25k_4541337', 'sx', 'sx56.tif'))

    path_network = gpd.read_file(os.path.join('Detailed-Path-Network', 'DARTMOOR NATIONAL PARK.gml'))

    tree = ET.parse(os.path.join('Detailed-Path-Network', 'DARTMOOR NATIONAL PARK.gml'))

    gpx_file = open(os.path.join('Walking routes', 'PFDartmoorwalk7BurratorReservoir.gpx'), 'r')
    g, path_nodes, path_network, nodes_id, node_coordinates = create_network(tree,path_network)
    bluebell_walk, coords, top_right, bottom_left = import_gpx(gpx_file)

    global_geom = []
    global_links = []
    global_nodes = []
    global_missing_links = []
    start_point = coords[0]
    for coord in coords[1:]:
        first_coordinate, first_node, last_coordinate, last_node = get_nearest_nodes(nodes_id, node_coordinates, start_point, coord, global_nodes, global_missing_links)
        path_gpd = gpx_to_path(first_node, last_node, g, path_network, global_geom, global_links)
        # plot_map(sheepstor_map, first_coordinate, last_coordinate, coords,
        #          bluebell_walk, path_network, path_nodes, path_gpd)
        start_point = coord
    global_path_gpd = gpd.GeoDataFrame({'fid': global_links, 'geometry': global_geom})
    global_nodes_gpd = gpd.GeoDataFrame({'geometry': global_nodes})

    global_missing_links_gpd = calculate_missing_sections(global_missing_links, global_path_gpd)
    plot_missing_links(global_missing_links_gpd, global_missing_links)


    plot_global_map(sheepstor_map, global_nodes_gpd, coords,
             bluebell_walk, path_network, path_nodes, global_path_gpd, top_right, bottom_left, global_missing_links_gpd)


if __name__ == '__main__':
    main()
