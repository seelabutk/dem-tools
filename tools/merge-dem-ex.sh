#!/bin/bash

gdal_merge.py -o merged.tiff chattanooga/chattanooga-e.DEM knoxville/knoxville-w.DEM corbin/corbin-e.DEM johnson_city/johnson_city-w.DEM
