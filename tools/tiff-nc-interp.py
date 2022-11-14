import os
import math
import argparse
import rasterio as rio
import numpy as np
import netCDF4 as nc
from tqdm import tqdm
from interp import *

# Parse command line arguments
parser = argparse.ArgumentParser(
    prog='tiff-nc-interp.py',
    description='Interpolate netCDF data to geoTiff file resolution, and crop geoTiff to bounds of netCDF data.'
)
parser.add_argument('-it', '--in-tiff', help='Path to geoTiff input file (Must have a single band)', required=True,
                    type=str)
parser.add_argument('-ot', '--out-tiff', help='Path to geoTiff output file.', required=True, type=str)
parser.add_argument('-in', '--in-nc', help='Path to netCDF input file.', required=True, type=str)
parser.add_argument('-on', '--out-nc', help='Path to netCDF output file.', required=True, type=str)
args = parser.parse_args()

############################################################################################################

# Open GeoTiff file
print(f'Reading input GeoTiff file [{args.in_tiff}]...')
dem = rio.open(args.in_tiff)
elevation = dem.read(1)
dem_lngs = np.array([dem.transform[2]], dtype=np.float64) + np.array(
    [i * np.array([dem.transform[0]], dtype=np.float64) for i in range(dem.width)], dtype=np.float64).flatten()
dem_lats = np.array([dem.transform[5]], dtype=np.float64) + np.array(
    [i * np.array([dem.transform[4]], dtype=np.float64) for i in range(dem.height)], dtype=np.float64).flatten()

# Open netCDF file
print(f'Reading input netCDF file [{args.in_nc}]...')
ds = nc.Dataset(args.in_nc)
nc_lats = ds['latitude'][:].astype(np.float64)
nc_lngs = ds['longitude'][:].astype(np.float64)

############################################################################################################

# Get Indices of GeoTiff data that are within the bounds of the netCDF data
print('Getting indices of GeoTiff data that are within the bounds of the netCDF data...')
nc_0_lat, nc_n_lat = nc_lats[0], nc_lats[-1]
nc_0_lng, nc_n_lng = nc_lngs[0], nc_lngs[-1]
dem_0_lat, dem_n_lat = dem_lats[0], dem_lats[-1]
dem_0_lng, dem_n_lng = dem_lngs[0], dem_lngs[-1]

# latitudes - decrease as i->n
min_lat_idx, max_lat_idx = -1, -1
for i in range(len(dem_lats)):
    lat = dem_lats[i]
    if lat < nc_0_lat and min_lat_idx == -1:  # correct
        min_lat_idx = i
    if lat > nc_n_lat:  # correct
        max_lat_idx = i

# longitudes - increase as i->n
min_lng_idx, max_lng_idx = -1, -1
for i in range(len(dem_lngs)):
    lng = dem_lngs[i]
    if lng < nc_n_lng:  # correct
        max_lng_idx = i
    if lng > nc_0_lng and min_lng_idx == -1:  # correct
        min_lng_idx = i

############################################################################################################
# create cropped geoTiff file

# crop the coords and elevation data
new_dem_lngs = dem_lngs[min_lng_idx:max_lng_idx + 1]
new_dem_lats = dem_lats[min_lat_idx:max_lat_idx + 1]
new_dem_elevation = elevation[min_lat_idx:max_lat_idx + 1, min_lng_idx:max_lng_idx + 1]

# build the geoTiff metadata
out_meta = dem.meta.copy()
out_transform = rio.transform.Affine(
    dem.meta['transform'][0], dem.meta['transform'][1], new_dem_lngs[0],
    dem.meta['transform'][3], dem.meta['transform'][4], new_dem_lats[0],
)
out_meta.update({'driver': 'GTiff',
                 'width': len(new_dem_lngs),
                 'height': len(new_dem_lats),
                 'dtype': 'int16',
                 'transform': out_transform})

# write file
print(f'Writing cropped GeoTiff file [{args.out_tiff}]...')
with rio.open(args.out_tiff, 'w', **out_meta) as dst:
    dst.write_band(1, new_dem_elevation.astype(rio.int16))


############################################################################################################

# min max scale
# def min_max_scale(data):
#     min_v, max_v = np.min(data), np.max(data)
#     if min_v == max_v:
#         return np.zeros(data.shape)
#     for i in range(len(data)):
#         for j in range(len(data[i])):
#             v = data[i][j]
#             data[i][j] = (v - min_v) / (max_v - min_v)
#     return data


############################################################################################################

# bilinear interpolation
# def blerp(q11, q12, q21, q22, x1, x2, y1, y2, x, y):
#     x2x1 = x2 - x1
#     y2y1 = y2 - y1
#     x2x = x2 - x
#     y2y = y2 - y
#     yy1 = y - y1
#     xx1 = x - x1
#     return 1 / (x2x1 * y2y1) * (
#             q11 * x2x * y2y +
#             q21 * xx1 * y2y +
#             q12 * x2x * yy1 +
#             q22 * xx1 * yy1
#     )
#
#
# # interpolate netCDF data for single timeslice for single variable to geoTiff resolution
# def interpolate(data):
#     data_new = np.zeros((max_lat_idx - min_lat_idx + 1, max_lng_idx - min_lng_idx + 1))
#     nc_lat_idx = 0
#     for dem_lat_idx in range(min_lat_idx, max_lat_idx + 1):
#         dem_lat = dem_lats[dem_lat_idx]
#         nc_lat, nc_lat_next = nc_lats[nc_lat_idx], nc_lats[nc_lat_idx + 1]
#         if dem_lat < nc_lat_next:
#             nc_lat_idx += 1
#             nc_lat, nc_lat_next = nc_lats[nc_lat_idx], nc_lats[nc_lat_idx + 1]
#         nc_lng_idx = 0
#         for dem_lng_idx in range(min_lng_idx, max_lng_idx + 1):
#             dem_lng = dem_lngs[dem_lng_idx]
#             nc_lng, nc_lng_next = nc_lngs[nc_lng_idx], nc_lngs[nc_lng_idx + 1]
#             if dem_lng > nc_lng_next:
#                 nc_lng_idx += 1
#                 nc_lng, nc_lng_next = nc_lngs[nc_lng_idx], nc_lngs[nc_lng_idx + 1]
#
#             data_new[dem_lat_idx - max_lat_idx, dem_lng_idx - max_lng_idx] = blerp(
#                 data[nc_lat_idx + 1][nc_lng_idx],  # Q11
#                 data[nc_lat_idx + 1][nc_lng_idx + 1],  # Q12
#                 data[nc_lat_idx][nc_lng_idx],  # Q21
#                 data[nc_lat_idx][nc_lng_idx + 1],  # Q22
#                 nc_lat_next,  # x1 - left x value
#                 nc_lat,  # x2 - right x value
#                 nc_lng,  # y1 - bottom y value
#                 nc_lng_next,  # y2 - top y value
#                 dem_lat,  # x - x interp value
#                 dem_lng,  # y - y interp value
#             )
#     return data_new


############################################################################################################

# Create interpolated netCDF file

# open dataset
print(f'Opening interpolated netCDF file [{args.out_nc}]...')
new_ds = nc.Dataset(args.out_nc, 'w', format="NETCDF4")

# create dimensions
new_ds_dim_lat = new_ds.createDimension("longitude", size=len(new_dem_lngs))
new_ds_dim_lon = new_ds.createDimension("latitude", size=len(new_dem_lats))
new_ds_dim_time = new_ds.createDimension("time", size=None)

# create lat lon variables
new_ds_var_lon = new_ds.createVariable('longitude', 'f4', ('longitude',))
new_ds_var_lat = new_ds.createVariable('latitude', 'f4', ('latitude',))
new_ds_var_lon[:] = new_dem_lngs
new_ds_var_lat[:] = new_dem_lats
from sklearn.preprocessing import MinMaxScaler
# begin interpolating data and writing it to the new netCDF file
print(f'Interpolating data and writing to [{args.out_nc}]...')
t_size = ds.dimensions['time'].size
i = 0
for k in tqdm(ds.variables.keys()):
    if k not in ['longitude', 'latitude', 'dem', 'bdy']:
        if k != 'Air_Pressure_at_Surface':
            continue
        new_ds_var = new_ds.createVariable(k, 'f4', ('time', 'latitude', 'longitude',))
        if 'units' in dir(ds[k]):
            new_ds_var.units = ds[k].units

        for t in tqdm(range(t_size), leave=False):
            data = ds[k][t,:].astype(np.float64)
            data = min_max_scale(data)
            data = interpolate(data.T,
                               max_lat_idx,
                               min_lat_idx,
                               max_lng_idx,
                               min_lng_idx,
                               dem_lats,
                               dem_lngs,
                               nc_lats,
                               nc_lngs)
            new_ds_var[t,:] = data
        i +=1
        if i > 2:
            break
# close the dataset
print(f'Closing interpolated netCDF file [{args.out_nc}]...')
new_ds.close()

print('Done!')
