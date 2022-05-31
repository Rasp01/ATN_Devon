from shapely.geometry import Point, mapping, LineString
import numpy as np
import rasterio
from rtree import index
import json
import networkx as nx
import geopandas as gpd
import os
import matplotlib.pyplot as plt



sheepstor_map = rasterio.open(os.path.join('OS Explorer Maps','Download_South+Dartmoor_2004150','raster-25k_4541337','sx','sx56.tif'))

print(sheepstor_map.crs)