#!/usr/bin/env python
# coding: utf-8

# Import packages
import os
import sys
import warnings
import matplotlib.pyplot as plt
#import numpy.ma as ma
import xarray as xr
import rioxarray as rxr
from shapely.geometry import mapping, box
import geopandas as gpd
import numpy as np
#import earthpy as et
#import earthpy.plot as ep

warnings.simplefilter('ignore')

# Open data with rioxarray
modis_lc = rxr.open_rasterio('/fmc_modis/land_cover.hdf') # MODIS veg cover data
modis_fmc = xr.open_dataset('/fmc_modis/modis_act_2020.nc')# MODIS FMC data
sf = gpd.read_file("C:/ACT_STATE_POLYGON_shp.shp") # ACT border shapefile

#Checking the coordinate systems
print(modis_lc.rio.crs)
print(modis_fmc.rio.crs)
print(sf.crs)

#Setting a coordinate system
modis_fmc=modis_fmc.rio.write_crs("EPSG:4283")
# Reprojecting the coordinates to make common coordinate system for all the dataset
modis_lc = modis_lc.rio.reproject("EPSG:3577", resolution=250.0) # reprojecting
modis_fmc = modis_fmc.rio.reproject("EPSG:3577", resolution=250.0)
sf = sf.to_crs("EPSG:3577")

# Cropping the image by ACT area
sf_geom = sf['geometry']
modis_lc_act = modis_lc.LC_Type1.rio.clip(sf_geom, from_disk=True)
modis_fmc_act = modis_fmc.fmc_mean.rio.clip(sf_geom, from_disk=True)

#Let's plot and see
f, ax = plt.subplots(1, 3, figsize=(25, 5))# let's plot both dataset
modis_lc_act.plot(ax=ax[0])
modis_fmc_act.sel(time = '2020-01-01').plot(ax=ax[1])
sf.plot(ax=ax[2], color="red")

np.set_printoptions(threshold=sys.maxsize) #we can plot everything in the numpy array
print(np.unique(modis_lc_act.values)) # printing unique values in the numpy array

#Creating a function to categorize veg cover, grass -1, shrub -2, forest -3
def categorize_lc(lc):
    return xr.where(lc == 1, 3, xr.where(lc == 2, 3, xr.where(lc == 3, 3, xr.where(lc == 4, 3, xr.where(lc == 5, 3, xr.where(lc ==6, 2,
        xr.where(lc ==7, 2, xr.where(lc ==8, 3, xr.where(lc ==9, 3, xr.where(lc ==10, 1, xr.where(lc ==12, 1, np.nan))))))))))) 

# applying a function
veg_mask = categorize_lc(modis_lc_act)
print(veg_mask)

#Converting data xr.array to xr.dataset
veg_mask = veg_mask.to_dataset()
modis_fmc_act = modis_fmc_act.to_dataset()

veg_mask = veg_mask.squeeze('band') # removing band from dimension
#veg_mask = veg_mask.expand_dims(time=modis.time) # we can also add time dimension if needed

time, y, x = modis_fmc_act.indexes.values() # Getting index values of fmc_modis
veg_mask_interp = veg_mask.interp(x=x, y=y, method="nearest") #interpolating x,y coordinates of veg_mask to make it overlap with fmc data
#This has been done because FMC and veg_mask didn't overlap perfectly, therefore I couldn't mask FMC based on veg_cover.

fmc_veg_cov = xr.merge([modis_fmc_act, veg_mask_interp])# merging fmc and veg_mask datset works perfectly after interpolating veg_mask x,y
#Exporting
fmc_veg_cov.to_netcdf('/fmc_veg_cov_2020.nc')

