from shapely.geometry import Point, LineString
import numpy as np
import rasterio
from rtree import index
import networkx as nx
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from cartopy import crs
import xml.etree.ElementTree as et
import gpxpy.gpx
from pyproj import CRS
from pyproj import Transformer


def create_network(tree):
    root = tree.getroot()
    startNodes = [startNode.attrib for startNode in
                  root.iter('{http://namespaces.os.uk/networks/detailedPathNetwork/1.0}startNode')]
    startNodes = [nodes['{http://www.w3.org/1999/xlink}href'] for nodes in startNodes]
    endNodes = [endNode.attrib for endNode in
                root.iter('{http://namespaces.os.uk/networks/detailedPathNetwork/1.0}endNode')]
    endNodes = [nodes['{http://www.w3.org/1999/xlink}href'] for nodes in endNodes]
    path_network = gpd.read_file(os.path.join('Detailed-Path-Network', 'DARTMOOR NATIONAL PARK.gml'))
    path_network['startNodes'] = startNodes
    path_network['endNodes'] = endNodes

    nodes_id = [nodes.attrib for nodes in
                root.iter('{http://namespaces.os.uk/networks/detailedPathNetwork/1.0}RouteNode')]
    nodes_id = [node['{http://www.opengis.net/gml/3.2}id'] for node in nodes_id]

    network_coords = [geom.text for geom in root.iter('{http://www.opengis.net/gml/3.2}pos')]
    node_coordinates = []
    for i in network_coords:
        coords = i.split(" ")
        coords = [float(j) for j in coords]
        node_coordinates.append(coords)

    node_geom = [Point(node) for node in node_coordinates]

    path_nodes = gpd.GeoDataFrame({'fid': nodes_id, 'geometry': node_geom})

    g = nx.Graph()
    for index, row in path_network.iterrows():
        g.add_edge(row['startNodes'], row['endNodes'], fid=row['gml_id'], weight=row['planimetricLength'])

    return g, path_nodes, path_network, nodes_id, node_coordinates


def import_gpx(gpx_file):
    gpx = gpxpy.parse(gpx_file)

    wsg84 = CRS.from_epsg(4326)
    osgb36 = CRS.from_epsg(27700)
    transformer = Transformer.from_crs(wsg84, osgb36)

    routes = gpx.routes
    waypoints = routes[0].points
    coords = []
    for waypoint in waypoints:
        coords.append(transformer.transform(waypoint.latitude, waypoint.longitude))
    bluebell_walk = gpd.GeoDataFrame(index=[0], crs='epsg:27700', geometry=[LineString(coords)])
    bluebell_walk.plot()
    return bluebell_walk, coords


def get_nearest_nodes(nodes_id, node_coordinates, start_point, coord):
    idx = index.Index()

    # set the bounds for the index
    for i in range(len(nodes_id)):
        left, bottom, right, top = (node_coordinates[i][0], node_coordinates[i][1],
                                    node_coordinates[i][0], node_coordinates[i][1])
        idx.insert(i, (left, bottom, right, top))

    for i in idx.nearest(start_point, 1):
        first_coordinate = node_coordinates[i]
        # need to add # so that it is identified in with the link
        first_node = ("#" + nodes_id[i])
    for i in idx.nearest(coord, 1):
        last_coordinate = node_coordinates[i]
        # need to add # so that it is identified in with the link
        last_node = ("#" + nodes_id[i])

    return first_coordinate, first_node, last_coordinate, last_node


def gpx_to_path(first_node, last_node, g, path_network):
    path = nx.dijkstra_path(g, first_node, last_node, weight='weight')

    geom = []
    links = []

    first_node = path[0]
    for node in path[1:]:
        link_fid = g.edges[first_node, node]['fid']
        links.append(link_fid)
        row = path_network.loc[path_network['gml_id'] == link_fid]
        geom.append(row['geometry'].cascaded_union)
        first_node = node

    path_gpd = gpd.GeoDataFrame({'fid': links, 'geometry': geom})
    path_gpd.plot()
    return path_gpd


def plot_map(sheepstor_map, first_coordinate, last_coordinate, coords,
             bluebell_walk, path_network, path_nodes, path_gpd):
    back_array = sheepstor_map.read(1)
    palette = np.array([value for key, value in sheepstor_map.colormap(1).items()])
    background_image = palette[back_array]
    bounds = sheepstor_map.bounds
    extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
    fig = plt.figure(figsize=(3, 3), dpi=500)
    ax = fig.add_subplot(1, 1, 1, projection=crs.OSGB())
    ax.imshow(background_image, origin='upper', extent=extent, zorder=0)
    plt.scatter(first_coordinate[0], first_coordinate[1], color='green', s=1, zorder=5)
    plt.scatter(last_coordinate[0], last_coordinate[1], color='green', s=1, zorder=5)
    bluebell_walk.plot(ax=ax, edgecolor='red', linewidth=0.5, zorder=2)
    plt.scatter(*zip(*coords), color='red', s=1, zorder=3)
    path_network.plot(ax=ax, edgecolor='blue', linewidth=0.5, zorder=4)
    path_nodes.plot(ax=ax, color='blue', markersize=1, zorder=4)
    path_gpd.plot(ax=ax, edgecolor='green', linewidth=0.5, zorder=5)
    display_extent = ((first_coordinate[0] - 500, first_coordinate[0] + 500,
                       first_coordinate[1] - 500, first_coordinate[1] + 500))
    ax.set_extent(display_extent, crs=crs.OSGB())
    plt.show()


def main():
    sheepstor_map = rasterio.open(
        os.path.join('OS Explorer Maps', 'Download_South+Dartmoor_2004150', 'raster-25k_4541337', 'sx', 'sx56.tif'))

    tree = et.parse(os.path.join('Detailed-Path-Network', 'DARTMOOR NATIONAL PARK.gml'))

    gpx_file = open(os.path.join('Walking routes', 'Bluebell walk.gpx'), 'r')
    g, path_nodes, path_network, nodes_id, node_coordinates = create_network(tree)
    bluebell_walk, coords = import_gpx(gpx_file)

    start_point = coords[0]
    for coord in coords[1:]:
        first_coordinate, first_node, last_coordinate, last_node = get_nearest_nodes(nodes_id, node_coordinates, start_point, coord)
        path_gpd = gpx_to_path(first_node, last_node, g, path_network)
        plot_map(sheepstor_map, first_coordinate, last_coordinate, coords,
                 bluebell_walk, path_network, path_nodes, path_gpd)
        start_point = coord


if __name__ == '__main__':
    main()
