import shapely.geometry
from shapely.ops import split, substring
from shapely.geometry import Point, mapping, LineString, Polygon, MultiPoint
import numpy as np
import rasterio
from rasterio import plot, mask
import networkx as nx
import pandas as pd
import geopandas as gpd
import os
import matplotlib.pyplot as plt
from cartopy import crs
from math import atan, degrees


def limit_elevation(elevation, study_area_shapely):
    study_area = mapping(study_area_shapely)
    elevation_mask, transform_index = mask.mask(elevation, [study_area], filled=False, crop=False)
    return elevation_mask, transform_index


def limit_dpn(path_network):
    OSDPN = path_network.intersection(study_area_shapely)
    OSDPN = OSDPN[~OSDPN.is_empty]
    osdpn_gpd = gpd.GeoDataFrame(geometry=OSDPN)
    return osdpn_gpd


def create_buffer(network_links, network_links_dpn, osdpn_gpd, buffer_value, buffer_cost):
    buffer = osdpn_gpd.buffer(buffer_value)
    buffer_gpd = gpd.GeoDataFrame(geometry=buffer)
    intersecting_links = gpd.sjoin(network_links_dpn, buffer_gpd)
    fid_list = intersecting_links.fid.tolist()
    network_links.loc[network_links['fid'].isin(fid_list), 'dpn'] = buffer_cost
    network_links_dpn.loc[network_links_dpn['fid'].isin(fid_list), 'dpn'] = buffer_cost
    network_links_dpn = network_links_dpn[network_links_dpn['dpn'] != buffer_cost]
    return network_links, network_links_dpn


def create_osdpn_buffers(network_links, osdpn_gpd):
    intersecting_nodes = gpd.sjoin(network_nodes, osdpn_gpd)
    intersecting_links = gpd.sjoin(network_links, osdpn_gpd)

    intersecting_links = intersecting_links.drop('index_right', 1)
    intersecting_links['dpn'] = 3
    dpn = intersecting_links["dpn"]
    network_links = network_links.join(dpn)
    network_links_dpn = network_links
    network_links_dpn = network_links_dpn[network_links_dpn['dpn'] != 3]



    network_links, network_links_dpn = create_buffer(network_links, network_links_dpn, osdpn_gpd, 40, 4)
    network_links, network_links_dpn = create_buffer(network_links, network_links_dpn, osdpn_gpd, 80, 5)
    network_links, network_links_dpn = create_buffer(network_links, network_links_dpn, osdpn_gpd, 120, 6)
    network_links, network_links_dpn = create_buffer(network_links, network_links_dpn, osdpn_gpd, 160, 7)

    return network_links


def create_intercept_points(osdpn_gpd):
    start_end_points = []
    for index1, row1 in osdpn_gpd.iterrows():
        if row1.geometry.geom_type != 'MultiLineString':
            path = row1.geometry
            # start_point = Point(path.coords[0])
            # end_point = Point(path.coords[-1])
            start_end_points.append(Point(path.coords[0]))
            start_end_points.append(Point(path.coords[-1]))

    start_end_points_gpd = gpd.GeoSeries(start_end_points)
    g = start_end_points_gpd.apply(lambda geom: geom.wkb)
    removed_duplicate_points = start_end_points_gpd.loc[g.drop_duplicates().index]
    osdpn_intercept_points = gpd.GeoDataFrame(geometry=removed_duplicate_points)
    osdpn_intercept_points['fid'] = range(1, len(osdpn_intercept_points) + 1)
    osdpn_intercept_points['fid'] = 'int_dpn' + osdpn_intercept_points['fid'].astype(str)
    return osdpn_intercept_points


def create_dpn_grid(osdpn_gpd, osdpn_intercept_points, network_links):
    geom_points = []
    geom_lines = []
    start_node = []
    end_node = []
    height = []
    length = []
    angles = []
    climb_time_forward = []
    dpn = []
    fid = []
    drop_fid = []
    count = 1
    for index1, row1 in osdpn_gpd.iterrows():
        if row1.geometry.geom_type != 'MultiLineString':
            line_coords = MultiPoint([Point(points) for points in row1.geometry.coords])
            result = split(row1.geometry, line_coords)
            previous_end_point = 0
            for linestring in result:

                intersects = {}
                local_height = []
                local_fid = []
                previous_intersecting_point = 0

                if previous_end_point != Point(linestring.coords[0]):

                    start_point = Point(linestring.coords[0])
                    intersects[count] = start_point
                    geom_points.append(start_point)
                    height.append(start_point.z)
                    local_height.append(start_point.z)

                    start_end = [Point(linestring.coords[0]), Point(linestring.coords[-1])]
                    start_end_gpd = gpd.GeoDataFrame(geometry=start_end)
                    point_join = gpd.sjoin(start_end_gpd, osdpn_intercept_points, how='left')

                    if pd.isna(point_join['fid'].values[0]):
                        start_point_fid = 'dpn_' + str(count)
                        fid.append(start_point_fid)
                        local_fid.append(start_point_fid)
                        count = count + 1
                    else:
                        start_point_fid = point_join['fid'].values[0]
                        fid.append(start_point_fid)
                        local_fid.append(start_point_fid)
                        count = count + 1

                    for index2, row2 in network_links.iterrows():
                        intersecting_point = linestring.intersection(row2.geometry)
                        if not intersecting_point.is_empty:
                            if intersecting_point != previous_intersecting_point:
                                fid_value = row2.fid
                                drop_fid.append(fid_value)
                                drop_fid.append(fid_value + 1)
                                network_node = [(Point(row2.geometry.coords[0])), (Point(row2.geometry.coords[-1]))]
                                if intersecting_point.geom_type != 'MultiPoint':

                                    intersects[count] = intersecting_point
                                    geom_points.append(intersecting_point)
                                    local_height.append(intersecting_point.z)
                                    height.append(intersecting_point.z)
                                    intersecting_point_fid = 'dpn_' + str(count)
                                    fid.append(intersecting_point_fid)
                                    local_fid.append(intersecting_point_fid)

                                    for node in network_node:

                                        # getting fid values for forward direction
                                        start_point = network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[
                                            0]
                                        end_point = intersecting_point_fid
                                        start_node.append(start_point)
                                        end_node.append(end_point)

                                        # get the reverse direction
                                        start_node.append(end_point)
                                        end_node.append(start_point)
                                        ## get line geometry
                                        line_forwards = LineString(
                                            [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                        line_backwards = LineString(
                                            [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                        geom_lines.append(line_forwards)
                                        geom_lines.append(line_backwards)

                                        # length
                                        length.append(line_forwards.length)
                                        length.append(line_backwards.length)

                                        # get angle of elevation change and climb time
                                        start_height_forward = end_height_backward = \
                                        network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                        # print(start_height)
                                        end_height_forward = start_height_backward = intersecting_point.z

                                        # print(end_height)
                                        if start_height_forward > end_height_forward:
                                            change_height = start_height_forward - end_height_forward
                                            climb_time_forward.append(0)
                                        else:
                                            change_height = end_height_forward - start_height_forward
                                            climb_time_forward.append(change_height / 10)

                                        # print(change_height)
                                        angle = degrees(atan((change_height / line_forwards.length)))
                                        # print(angle)
                                        angles.append(angle)
                                        dpn.append(2)

                                        if start_height_backward > end_height_backward:
                                            change_height = start_height_backward - end_height_backward
                                            climb_time_forward.append(0)
                                        else:
                                            change_height = end_height_backward - start_height_backward
                                            climb_time_forward.append(change_height / 10)

                                        # print(change_height)
                                        angle = degrees(atan((change_height / line_backwards.length)))
                                        # print(angle)
                                        angles.append(angle)
                                        dpn.append(2)

                                    count = count + 1
                                    previous_intersecting_point = intersecting_point
                                else:
                                    for pt in intersecting_point:
                                        intersects[count] = pt
                                        geom_points.append(pt)
                                        local_height.append(pt.z)
                                        height.append(pt.z)
                                        intersecting_point_fid = 'dpn_' + str(count)
                                        fid.append(intersecting_point_fid)
                                        local_fid.append(intersecting_point_fid)

                                        for node in network_node:

                                            # getting fid values for forward direction
                                            start_point = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[0]
                                            end_point = intersecting_point_fid
                                            start_node.append(start_point)
                                            end_node.append(end_point)

                                            # get the reverse direction
                                            start_node.append(end_point)
                                            end_node.append(start_point)
                                            ## get line geometry
                                            line_forwards = LineString(
                                                [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                            line_backwards = LineString(
                                                [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                            geom_lines.append(line_forwards)
                                            geom_lines.append(line_backwards)

                                            # length
                                            length.append(line_forwards.length)
                                            length.append(line_backwards.length)

                                            # get angle of elevation change and climb time
                                            start_height_forward = end_height_backward = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                            # print(start_height)
                                            end_height_forward = start_height_backward = intersecting_point.z

                                            # print(end_height)
                                            if start_height_forward > end_height_forward:
                                                change_height = start_height_forward - end_height_forward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_forward - start_height_forward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_forwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                            if start_height_backward > end_height_backward:
                                                change_height = start_height_backward - end_height_backward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_backward - start_height_backward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_backwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                        count = count + 1
                                        previous_intersecting_point = intersecting_point

                    end_point = Point(linestring.coords[-1])
                    intersects[count] = end_point
                    geom_points.append(end_point)
                    local_height.append(end_point.z)
                    height.append(end_point.z)

                    if pd.isna(point_join['fid'].values[1]):
                        start_point_fid = 'dpn_' + str(count)
                        fid.append(start_point_fid)
                        local_fid.append(start_point_fid)
                        count = count + 1
                    else:
                        start_point_fid = point_join['fid'].values[1]
                        fid.append(start_point_fid)
                        local_fid.append(start_point_fid)

                    previous_end_point = end_point
                    previous_end_point_fid = start_point_fid

                    points = MultiPoint(list(intersects.values()))
                    points_2d = shapely.wkb.loads(shapely.wkb.dumps(points, output_dimension=2))
                    line_2d = shapely.wkb.loads(shapely.wkb.dumps(linestring, output_dimension=2))
                    intersecting_points_2d = [Point(p.x, p.y) for p in points_2d]
                    intersects_gpd = gpd.GeoDataFrame(
                        {'fid': local_fid, 'height': local_height, 'geometry': intersecting_points_2d})
                    # adding the points to the linestring
                    line_coords = MultiPoint([Point(points) for points in line_2d.coords])
                    union = points_2d.union(line_coords)
                    # converting points to linestring
                    line = LineString(union)
                    # Splitting the line based on where intersections
                    splitted = split(line, points_2d)
                    for link in splitted:

                        # get start and end node
                        start_point = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'fid'].iloc[0]
                        end_point = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'fid'].iloc[0]
                        # start_node.append('dpn_' + str(start_point))
                        # end_node.append( 'dpn_' + str(end_point))
                        start_node.append(start_point)
                        end_node.append(end_point)

                        # get the reverse for backward direction
                        # start_node.append( 'dpn_' + str(end_point))
                        # end_node.append('dpn_' + str(start_point))
                        start_node.append(end_point)
                        end_node.append(start_point)

                        # get line geometry
                        geom_lines.append(link)
                        reverse_link = substring(link, link.length, 0)
                        geom_lines.append(reverse_link)

                        # length
                        length.append(link.length)
                        length.append(link.length)

                        # get angle of elevation change and climb time
                        start_height_forward = end_height_backward = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'height'].iloc[0]
                        # print(start_height)
                        end_height_forward = start_height_backward = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'height'].iloc[0]

                        # print(end_height)
                        if start_height_forward > end_height_forward:
                            change_height = start_height_forward - end_height_forward
                            climb_time_forward.append(0)
                        else:
                            change_height = end_height_forward - start_height_forward
                            climb_time_forward.append(change_height / 10)

                        # print(change_height)
                        angle = degrees(atan((change_height / link.length)))
                        # print(angle)
                        angles.append(angle)
                        dpn.append(0)

                        if start_height_backward > end_height_backward:
                            change_height = start_height_backward - end_height_backward
                            climb_time_forward.append(0)
                        else:
                            change_height = end_height_backward - start_height_backward
                            climb_time_forward.append(change_height / 10)

                        # print(change_height)
                        angle = degrees(atan((change_height / link.length)))
                        # print(angle)
                        angles.append(angle)
                        dpn.append(0)
                else:
                    start_point = previous_end_point
                    intersects[count] = start_point
                    local_height.append(start_point.z)
                    local_fid.append(previous_end_point_fid)
                    count = count + 1

                    start_end = [Point(linestring.coords[0]), Point(linestring.coords[-1])]
                    start_end_gpd = gpd.GeoDataFrame(geometry=start_end)
                    point_join = gpd.sjoin(start_end_gpd, osdpn_intercept_points, how='left')

                    for index2, row2 in network_links.iterrows():
                        intersecting_point = linestring.intersection(row2.geometry)
                        if not intersecting_point.is_empty:
                            if intersecting_point != previous_intersecting_point:
                                fid_value = row2.fid
                                drop_fid.append(fid_value)
                                drop_fid.append(fid_value + 1)
                                network_node = [(Point(row2.geometry.coords[0])), (Point(row2.geometry.coords[-1]))]
                                if intersecting_point.geom_type != 'MultiPoint':

                                    intersects[count] = intersecting_point
                                    geom_points.append(intersecting_point)
                                    local_height.append(intersecting_point.z)
                                    height.append(intersecting_point.z)
                                    intersecting_point_fid = 'dpn_' + str(count)
                                    fid.append(intersecting_point_fid)
                                    local_fid.append(intersecting_point_fid)

                                    for node in network_node:

                                        # getting fid values for forward direction
                                        start_point = network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[
                                            0]
                                        end_point = intersecting_point_fid
                                        start_node.append(start_point)
                                        end_node.append(end_point)

                                        # get the reverse direction
                                        start_node.append(end_point)
                                        end_node.append(start_point)
                                        ## get line geometry
                                        line_forwards = LineString(
                                            [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                        line_backwards = LineString(
                                            [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                        geom_lines.append(line_forwards)
                                        geom_lines.append(line_backwards)

                                        # length
                                        length.append(line_forwards.length)
                                        length.append(line_backwards.length)

                                        # get angle of elevation change and climb time
                                        start_height_forward = end_height_backward = \
                                        network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                        # print(start_height)
                                        end_height_forward = start_height_backward = intersecting_point.z

                                        # print(end_height)
                                        if start_height_forward > end_height_forward:
                                            change_height = start_height_forward - end_height_forward
                                            climb_time_forward.append(0)
                                        else:
                                            change_height = end_height_forward - start_height_forward
                                            climb_time_forward.append(change_height / 10)

                                        # print(change_height)
                                        angle = degrees(atan((change_height / line_forwards.length)))
                                        # print(angle)
                                        angles.append(angle)
                                        dpn.append(2)

                                        if start_height_backward > end_height_backward:
                                            change_height = start_height_backward - end_height_backward
                                            climb_time_forward.append(0)
                                        else:
                                            change_height = end_height_backward - start_height_backward
                                            climb_time_forward.append(change_height / 10)

                                        # print(change_height)
                                        angle = degrees(atan((change_height / line_backwards.length)))
                                        # print(angle)
                                        angles.append(angle)
                                        dpn.append(2)

                                    count = count + 1
                                    previous_intersecting_point = intersecting_point
                                else:
                                    for pt in intersecting_point:
                                        intersects[count] = pt
                                        geom_points.append(pt)
                                        local_height.append(pt.z)
                                        height.append(pt.z)
                                        intersecting_point_fid = 'dpn_' + str(count)
                                        fid.append(intersecting_point_fid)
                                        local_fid.append(intersecting_point_fid)

                                        for node in network_node:

                                            # getting fid values for forward direction
                                            start_point = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[0]
                                            end_point = intersecting_point_fid
                                            start_node.append(start_point)
                                            end_node.append(end_point)

                                            # get the reverse direction
                                            start_node.append(end_point)
                                            end_node.append(start_point)
                                            ## get line geometry
                                            line_forwards = LineString(
                                                [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                            line_backwards = LineString(
                                                [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                            geom_lines.append(line_forwards)
                                            geom_lines.append(line_backwards)

                                            # length
                                            length.append(line_forwards.length)
                                            length.append(line_backwards.length)

                                            # get angle of elevation change and climb time
                                            start_height_forward = end_height_backward = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                            # print(start_height)
                                            end_height_forward = start_height_backward = intersecting_point.z

                                            # print(end_height)
                                            if start_height_forward > end_height_forward:
                                                change_height = start_height_forward - end_height_forward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_forward - start_height_forward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_forwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                            if start_height_backward > end_height_backward:
                                                change_height = start_height_backward - end_height_backward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_backward - start_height_backward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_backwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                        count = count + 1
                                        previous_intersecting_point = intersecting_point

                    end_point = Point(linestring.coords[-1])
                    intersects[count] = end_point
                    geom_points.append(end_point)
                    local_height.append(end_point.z)
                    height.append(end_point.z)

                    if pd.isna(point_join['fid'].values[1]):
                        start_point_fid = 'dpn_' + str(count)
                        fid.append(start_point_fid)
                        local_fid.append(start_point_fid)
                        count = count + 1
                    else:
                        start_point_fid = point_join['fid'].values[1]
                        fid.append(start_point_fid)
                        local_fid.append(start_point_fid)

                    previous_end_point = end_point
                    previous_end_point_fid = start_point_fid

                    points = MultiPoint(list(intersects.values()))
                    points_2d = shapely.wkb.loads(shapely.wkb.dumps(points, output_dimension=2))
                    line_2d = shapely.wkb.loads(shapely.wkb.dumps(linestring, output_dimension=2))
                    intersecting_points_2d = [Point(p.x, p.y) for p in points_2d]
                    intersects_gpd = gpd.GeoDataFrame(
                        {'fid': local_fid, 'height': local_height, 'geometry': intersecting_points_2d})
                    # adding the points to the linestring
                    line_coords = MultiPoint([Point(points) for points in line_2d.coords])
                    union = points_2d.union(line_coords)
                    # converting points to linestring
                    line = LineString(union)
                    # Splitting the line based on where intersections
                    splitted = split(line, points_2d)
                    for link in splitted:

                        # get start and end node
                        start_point = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'fid'].iloc[0]
                        end_point = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'fid'].iloc[0]
                        # start_node.append('dpn_' + str(start_point))
                        # end_node.append( 'dpn_' + str(end_point))
                        start_node.append(start_point)
                        end_node.append(end_point)

                        # get the reverse for backward direction
                        # start_node.append( 'dpn_' + str(end_point))
                        # end_node.append('dpn_' + str(start_point))
                        start_node.append(end_point)
                        end_node.append(start_point)

                        # get line geometry
                        geom_lines.append(link)
                        reverse_link = substring(link, link.length, 0)
                        geom_lines.append(reverse_link)

                        # length
                        length.append(link.length)
                        length.append(link.length)

                        # get angle of elevation change and climb time
                        start_height_forward = end_height_backward = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'height'].iloc[0]
                        # print(start_height)
                        end_height_forward = start_height_backward = \
                        intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'height'].iloc[0]

                        # print(end_height)
                        if start_height_forward > end_height_forward:
                            change_height = start_height_forward - end_height_forward
                            climb_time_forward.append(0)
                        else:
                            change_height = end_height_forward - start_height_forward
                            climb_time_forward.append(change_height / 10)

                        # print(change_height)
                        angle = degrees(atan((change_height / link.length)))
                        # print(angle)
                        angles.append(angle)
                        dpn.append(0)

                        if start_height_backward > end_height_backward:
                            change_height = start_height_backward - end_height_backward
                            climb_time_forward.append(0)
                        else:
                            change_height = end_height_backward - start_height_backward
                            climb_time_forward.append(change_height / 10)

                        # print(change_height)
                        angle = degrees(atan((change_height / link.length)))
                        # print(angle)
                        angles.append(angle)
                        dpn.append(0)
        else:
            for single_linestring in row1.geometry:
                line_coords = MultiPoint([Point(points) for points in single_linestring.coords])
                result = split(single_linestring, line_coords)
                previous_end_point = 0
                for linestring in result:

                    intersects = {}
                    local_height = []
                    local_fid = []
                    previous_intersecting_point = 0

                    if previous_end_point != Point(linestring.coords[0]):

                        start_point = Point(linestring.coords[0])
                        intersects[count] = start_point
                        geom_points.append(start_point)
                        height.append(start_point.z)
                        local_height.append(start_point.z)

                        start_end = [Point(linestring.coords[0]), Point(linestring.coords[-1])]
                        start_end_gpd = gpd.GeoDataFrame(geometry=start_end)
                        point_join = gpd.sjoin(start_end_gpd, osdpn_intercept_points, how='left')

                        if pd.isna(point_join['fid'].values[0]):
                            start_point_fid = 'dpn_' + str(count)
                            fid.append(start_point_fid)
                            local_fid.append(start_point_fid)
                            count = count + 1
                        else:
                            start_point_fid = point_join['fid'].values[0]
                            fid.append(start_point_fid)
                            local_fid.append(start_point_fid)
                            count = count + 1

                        for index2, row2 in network_links.iterrows():
                            intersecting_point = linestring.intersection(row2.geometry)
                            if not intersecting_point.is_empty:
                                if intersecting_point != previous_intersecting_point:
                                    fid_value = row2.fid
                                    drop_fid.append(fid_value)
                                    drop_fid.append(fid_value + 1)
                                    network_node = [(Point(row2.geometry.coords[0])), (Point(row2.geometry.coords[-1]))]
                                    if intersecting_point.geom_type != 'MultiPoint':

                                        intersects[count] = intersecting_point
                                        geom_points.append(intersecting_point)
                                        local_height.append(intersecting_point.z)
                                        height.append(intersecting_point.z)
                                        intersecting_point_fid = 'dpn_' + str(count)
                                        fid.append(intersecting_point_fid)
                                        local_fid.append(intersecting_point_fid)

                                        for node in network_node:

                                            # getting fid values for forward direction
                                            start_point = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[0]
                                            end_point = intersecting_point_fid
                                            start_node.append(start_point)
                                            end_node.append(end_point)

                                            # get the reverse direction
                                            start_node.append(end_point)
                                            end_node.append(start_point)
                                            ## get line geometry
                                            line_forwards = LineString(
                                                [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                            line_backwards = LineString(
                                                [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                            geom_lines.append(line_forwards)
                                            geom_lines.append(line_backwards)

                                            # length
                                            length.append(line_forwards.length)
                                            length.append(line_backwards.length)

                                            # get angle of elevation change and climb time
                                            start_height_forward = end_height_backward = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                            # print(start_height)
                                            end_height_forward = start_height_backward = intersecting_point.z

                                            # print(end_height)
                                            if start_height_forward > end_height_forward:
                                                change_height = start_height_forward - end_height_forward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_forward - start_height_forward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_forwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                            if start_height_backward > end_height_backward:
                                                change_height = start_height_backward - end_height_backward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_backward - start_height_backward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_backwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                        count = count + 1
                                        previous_intersecting_point = intersecting_point
                                    else:
                                        for pt in intersecting_point:
                                            intersects[count] = pt
                                            geom_points.append(pt)
                                            local_height.append(pt.z)
                                            height.append(pt.z)
                                            intersecting_point_fid = 'dpn_' + str(count)
                                            fid.append(intersecting_point_fid)
                                            local_fid.append(intersecting_point_fid)

                                            for node in network_node:

                                                # getting fid values for forward direction
                                                start_point = \
                                                network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[0]
                                                end_point = intersecting_point_fid
                                                start_node.append(start_point)
                                                end_node.append(end_point)

                                                # get the reverse direction
                                                start_node.append(end_point)
                                                end_node.append(start_point)
                                                ## get line geometry
                                                line_forwards = LineString(
                                                    [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                                line_backwards = LineString(
                                                    [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                                geom_lines.append(line_forwards)
                                                geom_lines.append(line_backwards)

                                                # length
                                                length.append(line_forwards.length)
                                                length.append(line_backwards.length)

                                                # get angle of elevation change and climb time
                                                start_height_forward = end_height_backward = \
                                                network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                                # print(start_height)
                                                end_height_forward = start_height_backward = intersecting_point.z

                                                # print(end_height)
                                                if start_height_forward > end_height_forward:
                                                    change_height = start_height_forward - end_height_forward
                                                    climb_time_forward.append(0)
                                                else:
                                                    change_height = end_height_forward - start_height_forward
                                                    climb_time_forward.append(change_height / 10)

                                                # print(change_height)
                                                angle = degrees(atan((change_height / line_forwards.length)))
                                                # print(angle)
                                                angles.append(angle)
                                                dpn.append(2)

                                                if start_height_backward > end_height_backward:
                                                    change_height = start_height_backward - end_height_backward
                                                    climb_time_forward.append(0)
                                                else:
                                                    change_height = end_height_backward - start_height_backward
                                                    climb_time_forward.append(change_height / 10)

                                                # print(change_height)
                                                angle = degrees(atan((change_height / line_backwards.length)))
                                                # print(angle)
                                                angles.append(angle)
                                                dpn.append(2)

                                            count = count + 1
                                            previous_intersecting_point = intersecting_point

                        end_point = Point(linestring.coords[-1])
                        intersects[count] = end_point
                        geom_points.append(end_point)
                        local_height.append(end_point.z)
                        height.append(end_point.z)

                        if pd.isna(point_join['fid'].values[1]):
                            start_point_fid = 'dpn_' + str(count)
                            fid.append(start_point_fid)
                            local_fid.append(start_point_fid)
                            count = count + 1
                        else:
                            start_point_fid = point_join['fid'].values[1]
                            fid.append(start_point_fid)
                            local_fid.append(start_point_fid)

                        previous_end_point = end_point
                        previous_end_point_fid = start_point_fid

                        points = MultiPoint(list(intersects.values()))
                        points_2d = shapely.wkb.loads(shapely.wkb.dumps(points, output_dimension=2))
                        line_2d = shapely.wkb.loads(shapely.wkb.dumps(linestring, output_dimension=2))
                        intersecting_points_2d = [Point(p.x, p.y) for p in points_2d]
                        intersects_gpd = gpd.GeoDataFrame(
                            {'fid': local_fid, 'height': local_height, 'geometry': intersecting_points_2d})
                        # adding the points to the linestring
                        line_coords = MultiPoint([Point(points) for points in line_2d.coords])
                        union = points_2d.union(line_coords)
                        # converting points to linestring
                        line = LineString(union)
                        # Splitting the line based on where intersections
                        splitted = split(line, points_2d)
                        for link in splitted:

                            # get start and end node
                            start_point = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'fid'].iloc[0]
                            end_point = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'fid'].iloc[0]
                            # start_node.append('dpn_' + str(start_point))
                            # end_node.append( 'dpn_' + str(end_point))
                            start_node.append(start_point)
                            end_node.append(end_point)

                            # get the reverse for backward direction
                            # start_node.append( 'dpn_' + str(end_point))
                            # end_node.append('dpn_' + str(start_point))
                            start_node.append(end_point)
                            end_node.append(start_point)

                            # get line geometry
                            geom_lines.append(link)
                            reverse_link = substring(link, link.length, 0)
                            geom_lines.append(reverse_link)

                            # length
                            length.append(link.length)
                            length.append(link.length)

                            # get angle of elevation change and climb time
                            start_height_forward = end_height_backward = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'height'].iloc[0]
                            # print(start_height)
                            end_height_forward = start_height_backward = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'height'].iloc[0]

                            # print(end_height)
                            if start_height_forward > end_height_forward:
                                change_height = start_height_forward - end_height_forward
                                climb_time_forward.append(0)
                            else:
                                change_height = end_height_forward - start_height_forward
                                climb_time_forward.append(change_height / 10)

                            # print(change_height)
                            angle = degrees(atan((change_height / link.length)))
                            # print(angle)
                            angles.append(angle)
                            dpn.append(0)

                            if start_height_backward > end_height_backward:
                                change_height = start_height_backward - end_height_backward
                                climb_time_forward.append(0)
                            else:
                                change_height = end_height_backward - start_height_backward
                                climb_time_forward.append(change_height / 10)

                            # print(change_height)
                            angle = degrees(atan((change_height / link.length)))
                            # print(angle)
                            angles.append(angle)
                            dpn.append(0)
                    else:
                        start_point = previous_end_point
                        intersects[count] = start_point
                        local_height.append(start_point.z)
                        local_fid.append(previous_end_point_fid)
                        count = count + 1

                        start_end = [Point(linestring.coords[0]), Point(linestring.coords[-1])]
                        start_end_gpd = gpd.GeoDataFrame(geometry=start_end)
                        point_join = gpd.sjoin(start_end_gpd, osdpn_intercept_points, how='left')

                        for index2, row2 in network_links.iterrows():
                            intersecting_point = linestring.intersection(row2.geometry)
                            if not intersecting_point.is_empty:
                                if intersecting_point != previous_intersecting_point:
                                    fid_value = row2.fid
                                    drop_fid.append(fid_value)
                                    drop_fid.append(fid_value + 1)
                                    network_node = [(Point(row2.geometry.coords[0])), (Point(row2.geometry.coords[-1]))]
                                    if intersecting_point.geom_type != 'MultiPoint':

                                        intersects[count] = intersecting_point
                                        geom_points.append(intersecting_point)
                                        local_height.append(intersecting_point.z)
                                        height.append(intersecting_point.z)
                                        intersecting_point_fid = 'dpn_' + str(count)
                                        fid.append(intersecting_point_fid)
                                        local_fid.append(intersecting_point_fid)

                                        for node in network_node:

                                            # getting fid values for forward direction
                                            start_point = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[0]
                                            end_point = intersecting_point_fid
                                            start_node.append(start_point)
                                            end_node.append(end_point)

                                            # get the reverse direction
                                            start_node.append(end_point)
                                            end_node.append(start_point)
                                            ## get line geometry
                                            line_forwards = LineString(
                                                [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                            line_backwards = LineString(
                                                [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                            geom_lines.append(line_forwards)
                                            geom_lines.append(line_backwards)

                                            # length
                                            length.append(line_forwards.length)
                                            length.append(line_backwards.length)

                                            # get angle of elevation change and climb time
                                            start_height_forward = end_height_backward = \
                                            network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                            # print(start_height)
                                            end_height_forward = start_height_backward = intersecting_point.z

                                            # print(end_height)
                                            if start_height_forward > end_height_forward:
                                                change_height = start_height_forward - end_height_forward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_forward - start_height_forward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_forwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                            if start_height_backward > end_height_backward:
                                                change_height = start_height_backward - end_height_backward
                                                climb_time_forward.append(0)
                                            else:
                                                change_height = end_height_backward - start_height_backward
                                                climb_time_forward.append(change_height / 10)

                                            # print(change_height)
                                            angle = degrees(atan((change_height / line_backwards.length)))
                                            # print(angle)
                                            angles.append(angle)
                                            dpn.append(2)

                                        count = count + 1
                                        previous_intersecting_point = intersecting_point
                                    else:
                                        for pt in intersecting_point:
                                            intersects[count] = pt
                                            geom_points.append(pt)
                                            local_height.append(pt.z)
                                            height.append(pt.z)
                                            intersecting_point_fid = 'dpn_' + str(count)
                                            fid.append(intersecting_point_fid)
                                            local_fid.append(intersecting_point_fid)

                                            for node in network_node:

                                                # getting fid values for forward direction
                                                start_point = \
                                                network_nodes.loc[network_nodes['geometry'] == node, 'fid'].iloc[0]
                                                end_point = intersecting_point_fid
                                                start_node.append(start_point)
                                                end_node.append(end_point)

                                                # get the reverse direction
                                                start_node.append(end_point)
                                                end_node.append(start_point)
                                                ## get line geometry
                                                line_forwards = LineString(
                                                    [(node.x, node.y), (intersecting_point.x, intersecting_point.y)])
                                                line_backwards = LineString(
                                                    [(intersecting_point.x, intersecting_point.y), (node.x, node.y)])
                                                geom_lines.append(line_forwards)
                                                geom_lines.append(line_backwards)

                                                # length
                                                length.append(line_forwards.length)
                                                length.append(line_backwards.length)

                                                # get angle of elevation change and climb time
                                                start_height_forward = end_height_backward = \
                                                network_nodes.loc[network_nodes['geometry'] == node, 'height'].iloc[0]
                                                # print(start_height)
                                                end_height_forward = start_height_backward = intersecting_point.z

                                                # print(end_height)
                                                if start_height_forward > end_height_forward:
                                                    change_height = start_height_forward - end_height_forward
                                                    climb_time_forward.append(0)
                                                else:
                                                    change_height = end_height_forward - start_height_forward
                                                    climb_time_forward.append(change_height / 10)

                                                # print(change_height)
                                                angle = degrees(atan((change_height / line_forwards.length)))
                                                # print(angle)
                                                angles.append(angle)
                                                dpn.append(2)

                                                if start_height_backward > end_height_backward:
                                                    change_height = start_height_backward - end_height_backward
                                                    climb_time_forward.append(0)
                                                else:
                                                    change_height = end_height_backward - start_height_backward
                                                    climb_time_forward.append(change_height / 10)

                                                # print(change_height)
                                                angle = degrees(atan((change_height / line_backwards.length)))
                                                # print(angle)
                                                angles.append(angle)
                                                dpn.append(2)

                                            count = count + 1
                                            previous_intersecting_point = intersecting_point

                        end_point = Point(linestring.coords[-1])
                        intersects[count] = end_point
                        geom_points.append(end_point)
                        local_height.append(end_point.z)
                        height.append(end_point.z)

                        if pd.isna(point_join['fid'].values[1]):
                            start_point_fid = 'dpn_' + str(count)
                            fid.append(start_point_fid)
                            local_fid.append(start_point_fid)
                            count = count + 1
                        else:
                            start_point_fid = point_join['fid'].values[1]
                            fid.append(start_point_fid)
                            local_fid.append(start_point_fid)

                        previous_end_point = end_point
                        previous_end_point_fid = start_point_fid

                        points = MultiPoint(list(intersects.values()))
                        points_2d = shapely.wkb.loads(shapely.wkb.dumps(points, output_dimension=2))
                        line_2d = shapely.wkb.loads(shapely.wkb.dumps(linestring, output_dimension=2))
                        intersecting_points_2d = [Point(p.x, p.y) for p in points_2d]
                        intersects_gpd = gpd.GeoDataFrame(
                            {'fid': local_fid, 'height': local_height, 'geometry': intersecting_points_2d})
                        # adding the points to the linestring
                        line_coords = MultiPoint([Point(points) for points in line_2d.coords])
                        union = points_2d.union(line_coords)
                        # converting points to linestring
                        line = LineString(union)
                        # Splitting the line based on where intersections
                        splitted = split(line, points_2d)
                        for link in splitted:

                            # get start and end node
                            start_point = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'fid'].iloc[0]
                            end_point = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'fid'].iloc[0]
                            # start_node.append('dpn_' + str(start_point))
                            # end_node.append( 'dpn_' + str(end_point))
                            start_node.append(start_point)
                            end_node.append(end_point)

                            # get the reverse for backward direction
                            # start_node.append( 'dpn_' + str(end_point))
                            # end_node.append('dpn_' + str(start_point))
                            start_node.append(end_point)
                            end_node.append(start_point)

                            # get line geometry
                            geom_lines.append(link)
                            reverse_link = substring(link, link.length, 0)
                            geom_lines.append(reverse_link)

                            # length
                            length.append(link.length)
                            length.append(link.length)

                            # get angle of elevation change and climb time
                            start_height_forward = end_height_backward = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[0]), 'height'].iloc[0]
                            # print(start_height)
                            end_height_forward = start_height_backward = \
                            intersects_gpd.loc[intersects_gpd['geometry'] == Point(link.coords[-1]), 'height'].iloc[0]

                            # print(end_height)
                            if start_height_forward > end_height_forward:
                                change_height = start_height_forward - end_height_forward
                                climb_time_forward.append(0)
                            else:
                                change_height = end_height_forward - start_height_forward
                                climb_time_forward.append(change_height / 10)

                            # print(change_height)
                            angle = degrees(atan((change_height / link.length)))
                            # print(angle)
                            angles.append(angle)
                            dpn.append(0)

                            if start_height_backward > end_height_backward:
                                change_height = start_height_backward - end_height_backward
                                climb_time_forward.append(0)
                            else:
                                change_height = end_height_backward - start_height_backward
                                climb_time_forward.append(change_height / 10)

                            # print(change_height)
                            angle = degrees(atan((change_height / link.length)))
                            # print(angle)
                            angles.append(angle)
                            dpn.append(0)

    dpn_paths_nodes = gpd.GeoDataFrame({'fid': fid, 'height': height, 'dpn': 0, 'geometry': geom_points}, crs=27700)

    links_fid = range(1, len(geom_lines) + 1)
    dpn_paths_links = gpd.GeoDataFrame({'fid': links_fid, 'startnode': start_node,
                                        'endnode': end_node, 'length': length, 'angle': angles,
                                        'climb_time_forward': climb_time_forward, 'dpn': dpn, 'geometry': geom_lines},
                                       crs=27700)

    dpn_paths_links['fid'] = 'dpn_' + dpn_paths_links['fid'].astype(str)

    G = dpn_paths_nodes.geometry.apply(lambda geom: geom.wkb)
    dpn_paths_nodes = dpn_paths_nodes.loc[G.drop_duplicates().index]

    G = dpn_paths_links.geometry.apply(lambda geom: geom.wkb)
    dpn_paths_links = dpn_paths_links.loc[G.drop_duplicates().index]

    return dpn_paths_nodes, dpn_paths_links, drop_fid


def intergrate_dpn(dpn_paths_nodes, network_nodes, dpn_paths_links, network_links, drop_fid):
    network_links = network_links[~network_links.fid.isin(drop_fid)]
    network_links['fid'] = 'al_' + network_links['fid'].astype(str)
    network_nodes = gpd.GeoDataFrame(pd.concat([network_nodes, dpn_paths_nodes], ignore_index=True))
    network_links = gpd.GeoDataFrame(pd.concat([network_links, dpn_paths_links], ignore_index=True))

    network_links["dpn"] = network_links["dpn"].fillna(8)
    network_nodes["dpn"] = network_nodes["dpn"].fillna(8)
    return network_nodes, network_links


def access_land_removal(access_land, network_nodes, network_links):
    access_land = access_land.to_crs('EPSG:27700')
    access_land_intersection = access_land.intersection(study_area_shapely)
    access_land_intersection = access_land_intersection[~access_land_intersection.is_empty]
    access_land_gpd = gpd.GeoDataFrame(geometry=access_land_intersection)

    intersecting_nodes = gpd.sjoin(network_nodes, access_land_gpd)
    intersecting_links = gpd.sjoin(network_links, access_land_gpd)

    outside_nodes = network_nodes.drop(intersecting_nodes.index)
    outside_links = network_links.drop(intersecting_links.index)

    intersecting_nodes_removed = outside_nodes.query('dpn != 0')
    intersecting_links_removed = outside_links.query('dpn != 0')

    network_nodes = network_nodes.drop(intersecting_nodes_removed.index)
    network_links = network_links.drop(intersecting_links_removed.index)
    return network_nodes, network_links


def remove_from_network(all_obstructions, network_nodes, network_links):
    obstructions = all_obstructions.intersection(study_area_shapely)
    obstructions = obstructions[~obstructions.is_empty]
    obstructions_gpd = gpd.GeoDataFrame(geometry=obstructions)

    intersecting_nodes = gpd.sjoin(network_nodes, obstructions_gpd)
    intersecting_links = gpd.sjoin(network_links, obstructions_gpd)
    intersecting_links = intersecting_links.drop(intersecting_links[(intersecting_links.dpn == 0)].index)
    intersecting_nodes = intersecting_nodes.drop(intersecting_nodes[(intersecting_nodes.dpn == 0)].index)

    network_nodes = network_nodes[~ network_nodes.isin(intersecting_nodes)].dropna()
    cond = network_links['fid'].isin(intersecting_links['fid'])
    network_links.drop(network_links[cond].index, inplace=True)
    return network_nodes, network_links


def master_map_removal(master_map_polygons, master_map_lines, network_nodes, network_links):
    # polygon
    all_obstructions = master_map_polygons.query(
        'descriptiveterm =="Slope" or descriptiveterm =="Cliff" or descriptivegroup =="Inland Water"')
    network_nodes, network_links = remove_from_network(all_obstructions, network_nodes, network_links)
    # lines
    all_obstructions = master_map_lines.query(
        'physicalpresence =="Obstructing" or descriptiveterm =="Bottom Of Slope" or descriptiveterm =="Bottom Of Cliff" '
        'or descriptiveterm =="Top Of Slope" or descriptiveterm =="Top Of Cliff"or descriptiveterm =="Watercourse"')
    network_nodes, network_links = remove_from_network(all_obstructions, network_nodes, network_links)

    return network_nodes, network_links


def drop_dangerous_slope(network_links):
    network_links = network_links.drop(network_links[(network_links.angle >= 30)].index)
    return network_links


def land_cover_classification(land_use, network_nodes, network_links, study_area):
    land_use_mask, land_use_transform_index = mask.mask(land_use, [study_area], filled=False, crop=False)

    # for links
    land_use_value = []
    for index, row in network_links.iterrows():
        location = (row.geometry.centroid.x, row.geometry.centroid.y)
        row, col = land_use.index(location[0], location[1])
        # print("Point Corresponds to row, col: %d, %d"%(row,col))
        # print("Raster value on point %.2f \n"%land_use.read(1)[row,col])
        land_use_value.append(land_use.read(1)[row, col])
    network_links['land_use'] = land_use_value

    network_links = network_links.drop(network_links[(network_links.land_use == 29)].index)

    surface_cost = []
    for index, row in network_links.iterrows():
        if row.land_use == 1 or row.land_use == 3 or row.land_use == 9 or row.land_use == 26 or row.land_use == 27:
            surface_cost.append(0)
        elif row.land_use == 6 or row.land_use == 7 or row.land_use == 2:
            surface_cost.append(1)
        elif row.land_use == 10 or row.land_use == 11 or row.land_use == 12 or row.land_use == 303 \
                or row.land_use == 314 or row.land_use == 328:
            surface_cost.append(2)
        elif row.land_use == 13 or row.land_use == 14 or row.land_use == 19 or row.land_use == 23:
            surface_cost.append(3)
        else:
            surface_cost.append(4)
    network_links["surface_cost"] = surface_cost

    # for nodes
    land_use_value = []
    for index, row in network_nodes.iterrows():
        location = (row.geometry.x, row.geometry.y)
        row, col = land_use.index(location[0], location[1])
        # print("Point Corresponds to row, col: %d, %d"%(row,col))
        # print("Raster value on point %.2f \n"%land_use.read(1)[row,col])
        land_use_value.append(land_use.read(1)[row, col])
    network_nodes['land_use'] = land_use_value

    return network_nodes, network_links


def naismiths_rule(network_links):
    total_time = []
    for index, row in network_links.iterrows():
        length_time = row.length / (500 / 6)
        total_time.append(length_time + row.climb_time_forward)

    network_links['total_time'] = total_time
    return network_links


def create_graph(network_links):
    graph = nx.DiGraph()
    for index, row in network_links.iterrows():
        graph.add_edge(row['startnode'], row['endnode'], fid=row['fid'], length=row.length, time=row.total_time)
    return graph


def create_paths(graph):
    # points = ["al_6", "al_1670"]
    points = ["int_dpn27", "dpn_271"]

    # get the shortest path with time weight
    weighted_path_forward = nx.dijkstra_path(graph, source=points[0], target=points[1], weight='time')
    weighted_path_backward = nx.dijkstra_path(graph, source=points[1], target=points[0], weight='time')
    return weighted_path_forward, weighted_path_backward


def create_path_gpd(graph, weighted_path, network_links):
    geom = []
    links = []
    first_node = weighted_path[0]
    for node in weighted_path[1:]:
        link_fid = graph.edges[first_node, node]['fid']
        links.append(link_fid)
        row = network_links.loc[network_links['fid'] == link_fid]
        geom.append(row['geometry'].cascaded_union)
        first_node = node

    weighted_path_gpd = gpd.GeoDataFrame({'fid': links, 'geometry': geom})
    return weighted_path_gpd


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
study_area = mapping(study_area_shapely)

SX77_map = rasterio.open(
    os.path.join('OS Explorer Maps', 'Download_SX77-Haytor_2033809', 'raster-25k_4596071', 'sx', 'sx77.tif'))

# elevation = rasterio.open(os.path.join('OS Elevation','SX77_elevation','terrain-5-dtm_4616587','sx','sx77ne_nw','w001001.adf'))

# network_links = gpd.read_file(os.path.join('Study_area', 'SX7677', 'Final Networks', 'network_links_al.geojson'))
# network_nodes = gpd.read_file(os.path.join('Study_area', 'SX7677', 'Final Networks', 'network_nodes_al.geojson'))

network_links = gpd.read_file(os.path.join('Study_area', 'SX7478', 'network_links_al.geojson'))
network_nodes = gpd.read_file(os.path.join('Study_area', 'SX7478', 'network_nodes_al.geojson'))

path_network = gpd.read_file(os.path.join('Detailed-Path-Network', 'DARTMOOR NATIONAL PARK.gml'), layer='RouteLink')

access_land = gpd.read_file(
    os.path.join('Devon County Council', 'Public_Access', 'Access_Land_and_Dartmoor_Commons.geojson'))

master_map_polygons = gpd.read_file(
    os.path.join('MasterMap', 'Download_SX77+-+Haytor_2035912', 'mastermap-topo_4600023',
                 'mastermap-topo_4600023.gpkg'), layer="Topographicarea")

master_map_lines = gpd.read_file(os.path.join('MasterMap', 'Download_SX77+-+Haytor_2035912', 'mastermap-topo_4600023',
                                              'mastermap-topo_4600023.gpkg'), layer="Topographicline")

land_use = rasterio.open(os.path.join('Environmental Land Use Data', 'Dartmoor_land_use_converted.tif'))

# elevation_mask,transform_index = limit_elevation(elevation, study_area_shapely)


osdpn_gpd = limit_dpn(path_network)

network_links = create_osdpn_buffers(network_links, osdpn_gpd)

osdpn_intercept_points = create_intercept_points(osdpn_gpd)

dpn_paths_nodes, dpn_paths_links, drop_fid = create_dpn_grid(osdpn_gpd, osdpn_intercept_points, network_links)

# dpn_paths_nodes.to_file("Study_area/SX7677/Final Networks/dpn_paths_nodes.geojson", driver='GeoJSON', crs='EPSG:27700')
# dpn_paths_links.to_file("Study_area/SX7677/Final Networks/dpn_paths_links.geojson", driver='GeoJSON', crs='EPSG:27700')

# dpn_paths_nodes.to_file("Study_area/SX7478/dpn_paths_nodes.geojson", driver='GeoJSON', crs='EPSG:27700')
# dpn_paths_links.to_file("Study_area/SX7478/dpn_paths_links.geojson", driver='GeoJSON', crs='EPSG:27700')

network_nodes, network_links = intergrate_dpn(dpn_paths_nodes, network_nodes, dpn_paths_links, network_links, drop_fid)

network_nodes, network_links = access_land_removal(access_land, network_nodes, network_links)

network_nodes, network_links = master_map_removal(master_map_polygons, master_map_lines, network_nodes, network_links)

network_links = drop_dangerous_slope(network_links)

network_nodes, network_links = land_cover_classification(land_use, network_nodes, network_links, study_area)

network_links = naismiths_rule(network_links)

# network_links.to_file("Study_area/SX7677/Final Networks/network_links_dpn.geojson", driver='GeoJSON', crs='EPSG:27700')
# network_nodes.to_file("Study_area/SX7677/Final Networks/network_nodes_dpn.geojson", driver='GeoJSON', crs='EPSG:27700')

# network_links.to_file("Study_area/SX7478/network_links_dpn.geojson", driver='GeoJSON', crs='EPSG:27700')
# network_nodes.to_file("Study_area/SX7478/network_nodes_dpn.geojson", driver='GeoJSON', crs='EPSG:27700')

graph = create_graph(network_links)

weighted_path_forward, weighted_path_backward = create_paths(graph)

weighted_path_forward_gpd = create_path_gpd(graph, weighted_path_forward, network_links)
weighted_path_backward_gpd = create_path_gpd(graph, weighted_path_backward, network_links)

plot_network(SX77_map, study_area_shapely
             , network_nodes, network_links, weighted_path_forward_gpd, weighted_path_backward_gpd)
