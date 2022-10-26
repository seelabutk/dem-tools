
import sys
from osgeo import gdal

import netCDF4 as nc
import numpy as np


in_tiff = gdal.Open(sys.argv[1])

nc_file = nc.Dataset(sys.argv[2])

delta = 0.001

max_lat = np.max(nc_file['latitude'][:])
min_lat = np.min(nc_file['latitude'][:])
max_lon = np.max(nc_file['longitude'][:])
min_lon = np.min(nc_file['longitude'][:])

print(f'lat: ({min_lat}, {max_lat})')
print(f'lon: ({min_lon}, {max_lon})')

# path to where you want the clipped raster
out_tiff = sys.argv[3]

ds = gdal.Translate(out_tiff , in_tiff,
                    projWin = [min_lon-delta, max_lat+delta,
                               max_lon+delta, min_lat-delta]) # OR [ulx, uly, lrx, lry]
ds = None