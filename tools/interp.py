import numpy as np
from numba import njit

# min max scale

# @profile
@njit(fastmath=True)
def min_max_scale(data):
    min_v, max_v = np.min(data), np.max(data)
    if min_v == max_v:
        return np.zeros(data.shape)
    for i in range(len(data)):
        for j in range(len(data[i])):
            v = data[i][j]
            data[i][j] = (v - min_v) / (max_v - min_v)
    return data

# @profile
@njit(fastmath=True)
def blerp(q11,q12,q21,q22,x1,x2,y1,y2,x,y):
    x2x1 = x2 - x1
    y2y1 = y2 - y1
    x2x = x2 - x
    y2y = y2 - y
    yy1 = y - y1
    xx1 = x - x1
    return 1 / (x2x1*  y2y1) * (
            q11 * x2x * y2y +
            q21 * xx1 * y2y +
            q12 * x2x * yy1 +
            q22 * xx1 * yy1
    )

# @profile
@njit(fastmath=True)
def interpolate(data,
                max_lat_idx,
                min_lat_idx,
                max_lng_idx,
                min_lng_idx,
                dem_lats,
                dem_lngs,
                nc_lats,
                nc_lngs):
    data_new = np.zeros((max_lat_idx - min_lat_idx + 1, max_lng_idx - min_lng_idx + 1))
    nc_lat_idx = 0
    for dem_lat_idx in range(min_lat_idx,max_lat_idx+1):
        # for dem_lat_idx in range(min_lat_idx,max_lat_idx+1):
        dem_lat = dem_lats[dem_lat_idx]
        nc_lat, nc_lat_next = nc_lats[nc_lat_idx], nc_lats[nc_lat_idx+1]
        if dem_lat < nc_lat_next:
            nc_lat_idx += 1
            nc_lat, nc_lat_next = nc_lats[nc_lat_idx], nc_lats[nc_lat_idx+1]
        nc_lng_idx = 0
        for dem_lng_idx in range(min_lng_idx, max_lng_idx+1):
            dem_lng = dem_lngs[dem_lng_idx]
            nc_lng, nc_lng_next = nc_lngs[nc_lng_idx], nc_lngs[nc_lng_idx+1]
            if dem_lng > nc_lng_next:
                nc_lng_idx += 1
                nc_lng, nc_lng_next = nc_lngs[nc_lng_idx], nc_lngs[nc_lng_idx+1]
            data_new[dem_lat_idx-max_lat_idx, dem_lng_idx-max_lng_idx] = blerp(
                data[nc_lat_idx+1][nc_lng_idx], # Q11 - bottom right value
                data[nc_lat_idx+1][nc_lng_idx+1], # Q12 - top right value
                data[nc_lat_idx][nc_lng_idx], # 21 - bottom left value
                data[nc_lat_idx][nc_lng_idx+1], # Q22 - top left value Q

                nc_lat_next, # x1 - left x value
                nc_lat, # x2 - right x value

                nc_lng, # y1 - bottom y value
                nc_lng_next, # y2 - top y value

                dem_lat, # x - x interp value
                dem_lng,  # y - y interp value
            )
    return data_new


