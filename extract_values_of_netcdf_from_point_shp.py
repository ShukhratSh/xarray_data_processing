
import numpy as np
import xarray as xr
import rioxarray
import geopandas as gpd

# Plotting options
sns.set(font_scale=1.3)
sns.set_style("white")

# reading MODIS FMC data between 2014 and 2020
fmc_modis_2014 = xr.open_dataset("D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2014.nc")
fmc_modis_2015 = xr.open_dataset("D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2015.nc")
fmc_modis_2016 = xr.open_dataset("D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2016.nc")
fmc_modis_2017 = xr.open_dataset("D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2017.nc")
fmc_modis_2018 = xr.open_dataset("D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2018.nc")
fmc_modis_2019 = xr.open_dataset("D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2019.nc")
fmc_modis_2020 = xr.open_dataset("D:/fmc_modis/fmc_modis_lc/fmc_veg_cov_2020.nc")

# Concatenating FMC data from different times
fmc_lut_modis = xr.concat([fmc_modis_2014, fmc_modis_2015, fmc_modis_2016, fmc_modis_2017, fmc_modis_2018, fmc_modis_2019, fmc_modis_2020], dim='time')

# Slicing the data between 2015-10-01 - 2020-12-31
fmc_lut_modis = fmc_lut_modis.sel(time=slice('2015-10-01','2020-12-31')) # selecting the data between 2015-2020 to match with Sentinel data

#Reprojecting and resampling
fmc_lut_modis = fmc_lut_modis.rio.reproject("EPSG:3577", resolution=500.0) # resampling MODIS data to 500 m spatial resolution 

#Reading a point shapefile
honeysuckle_shp = gpd.read_file("/honeysuckle.shp")

#Extracting x, y coordinate values
x = honeysuckle_shp.geometry.x
y = honeysuckle_shp.geometry.y

# Extracting the time-series FMC data for the point location of Honeysuckle
ds_honeysuckle = fmc_lut_modis.sel(x=x, y=y, method="nearest")
# Save the data
ds_honeysuckle.to_netcdf("/honeysuckle.nc")
