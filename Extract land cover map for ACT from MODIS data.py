#!/usr/bin/env python
# coding: utf-8

# In[159]:


# Import packages
import os
import warnings
import matplotlib.pyplot as plt
import numpy.ma as ma
import xarray as xr
import rioxarray as rxr
from shapely.geometry import mapping, box
import geopandas as gpd
import earthpy as et
import earthpy.plot as ep
import matplotlib.pyplot as plt

warnings.simplefilter('ignore')


# In[160]:


# Open data with rioxarray
modis_lc = rxr.open_rasterio('D:/fmc_modis/land_cover/MCD12Q1.A2019001.h30v12.006.2020212133248.hdf')
modis_fmc = xr.open_dataset('D:/fmc_modis/modis_act_2020.nc')
sf = gpd.read_file("C:/Users/u6262380/Downloads/act_state_polygon_shp/ACT_STATE_POLYGON_shp.shp")


# In[161]:


#Checking the coordinate systems
print(modis_lc.rio.crs)
print(modis_fmc.rio.crs)
print(sf.crs)


# In[162]:


#Setting a coordinate system
modis_fmc=modis_fmc.rio.write_crs("EPSG:4283")


# In[163]:


print(modis_fmc.rio.crs)
print(modis_lc.rio.crs)
print(sf.crs)


# In[164]:


modis_lc = modis_lc.rio.reproject("EPSG:3577", resolution=250.0) # reprojecting
modis_fmc = modis_fmc.rio.reproject("EPSG:3577", resolution=250.0)
sf = sf.to_crs("EPSG:3577")


# In[165]:


# Cropping the image by ACT area
sf_geom = sf['geometry']
modis_lc_act = modis_lc.LC_Type1.rio.clip(sf_geom, from_disk=True)
modis_fmc_act = modis_fmc.fmc_mean.rio.clip(sf_geom, from_disk=True)


# In[166]:


#Let's plot and see
f, ax = plt.subplots(1, 3, figsize=(25, 5))# let's plot both dataset
modis_lc_act.plot(ax=ax[0])
modis_fmc_act.sel(time = '2020-01-01').plot(ax=ax[1])
sf.plot(ax=ax[2], color="red")


# In[167]:


#modis_lc_act.to_netcdf('D:/fmc_modis/lc_modis_1.tif')
#modis_fmc_act.sel(time = '2014-01-01').to_netcdf('D:/fmc_modis/fmc_modis_1.tif')


# In[168]:


import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
print(np.unique(modis_lc_act.values))


# In[169]:


#Creating a function to categorize veg cover, grass -1, shrub -2, forest -3
def categorize_lc(lc):
    return xr.where(lc == 1, 3, xr.where(lc == 2, 3, xr.where(lc == 3, 3, xr.where(lc == 4, 3, xr.where(lc == 5, 3, xr.where(lc ==6, 2,
        xr.where(lc ==7, 2, xr.where(lc ==8, 3, xr.where(lc ==9, 3, xr.where(lc ==10, 1, xr.where(lc ==12, 1, np.nan))))))))))) 


# In[170]:


# applying a function
veg_mask = categorize_lc(modis_lc_act)
print(veg_mask)


# In[171]:


#Converting data xr.array to xr.dataset
veg_mask = veg_mask.to_dataset()
modis_fmc_act = modis_fmc_act.to_dataset()


# In[172]:


veg_mask = veg_mask.squeeze('band') # removing band from dimension
#veg_mask = veg_mask.expand_dims(time=modis.time) # adding time dimension
veg_mask


# In[173]:


time, y, x = modis_fmc_act.indexes.values() # Getting index values of fmc_modis
veg_mask_interp = veg_mask.interp(x=x, y=y, method="nearest") #interpolating x,y coordinates of veg_mask to make it overlap with fmc data
#This has been done because FMC aand veg_mask derived from MOdis didn't overlap perfectly, therefore I couldn't mask FMC based on veg_cover.


# In[175]:


fmc_veg_cov = xr.merge([modis_fmc_act, veg_mask_interp])# merging fmc and coord. interpolated veg_mask datset
fmc_veg_cov


# In[176]:


fmc_veg_cov.to_netcdf('D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2020.nc')

