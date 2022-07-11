"""

Author: Marine Remaud
Reads the forest data from the Pucher inventory
Finds the latitudes and longitudes coordinates (wgs 84 system) from the original coordinate espg
Converts the data from the tif format to a pandas dataframe, than creates a netcdf

"""

from osgeo import gdal,osr,ogr
import pyproj
import xarray as xr
from pyproj import Proj, transform
pyproj.datadir.set_data_dir("/opt/anaconda2/envs/test_esmpy/share/proj")
import geopandas as gpd
import rasterio
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#-----------------READ THE DIAMETER----------------------------------------------------------
parameter="Diameter"
homedata="/home/surface1/mremaud/FOREST"                       #Directory of the ForestInventory (Pucher data)
path=homedata+"/ForestInventory/"+parameter+"/dbh.tif"         #Mean diameter over the 8x8 km grid
path_sd=homedata+"/ForestInventory/"+parameter+"/dbh_sd.tif"   #Standard deviation of the diameter over the grid 8x8 km 

#Read the coordinate espg (not used after)
d=gdal.Open(path) #Open the file dbh.tif with gdal
proj=osr.SpatialReference(wkt=d.GetProjection())
espg_in=proj.GetAttrValue("AUTHORITY",1) #Find the reference of the coordinate system

# get the existing coordinate system
old_cs= osr.SpatialReference()
old_cs.ImportFromWkt(d.GetProjectionRef())

# create the new coordinate system
wgs84_wkt = """
GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]"""
new_cs = osr.SpatialReference()
new_cs.ImportFromWkt(wgs84_wkt)

# create a transform object to convert between coordinate systems
transform = osr.CoordinateTransformation(old_cs,new_cs) 

#Create two netcdf files (for the mean diameter and its standard deviation) with the old coordinate system
#The file contains only one variable of dimension number of longitudes x number of latitudes
#We create a netcdf here to get the x and y coordinates 
gdal.Translate(parameter+".nc",path,format="NetCDF")
gdal.Translate(parameter+"_sd.nc",path_sd,format="NetCDF")

#Open the netcdf
Val=xr.open_dataset(parameter+".nc")
Val_sd=xr.open_dataset(parameter+"_sd.nc")

#Determine the longitudes (coo_lon) and latitudes (coo_lon) of each observed location
x1,y1 = Val.x.values,Val.y.values
coo_lat=np.zeros((len(y1),len(x1)))
coo_lon=np.zeros((len(y1),len(x1)))

for ix,xx in enumerate(x1):
 for iy,yy in enumerate(y1): 
   latlong = transform.TransformPoint(xx,yy)
   coo_lat[iy,ix]=latlong[0]
   coo_lon[iy,ix]=latlong[1]

#Contour map of the mean diameter---------------------------------------------------------------------------
latmin=35
latmax=60
lonmin=-9
lonmax=29

m = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,\
              llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='c')
m.drawcoastlines()
n_graticules = 18
parallels = np.arange(-90., 90, 30)
meridians = np.arange(0., 360., 60)
lw = 0.5
dashes = [1,1] # 5 dots, 7 spaces... repeat
graticules_color = 'grey'
m.drawparallels(parallels, linewidth=0.5,labels=[1,0,0,0], dashes=dashes, color=graticules_color,zorder=20)
m.drawmeridians(meridians,linewidth=0.5, labels=[0,0,0,1],dashes=dashes, color=graticules_color,zorder=20)
Y, X = m(coo_lon,coo_lat)
ax=m.contourf(Y,X,Val.Band1.values,cmap=plt.cm.PiYG_r,extend='both')
#---------------------------------------------------------------------------------------------------------

#Create new dataframe

#Fixed and Indicative longitude and latitude (they are not correct)
#The correct longitudes and latitudes are contained in the variables coo_lat, coo_lon
vecx=[]
vecy=[]
for xx in x1:
  yy=y1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecx.append(latlong[1])  
for yy in y1:
  xx=x1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecy.append(latlong[0])
vecy=np.round(np.asarray(vecy),2)
vecx=np.round(np.asarray(vecx),2)


#Store the variable in xarray (Val.Band1.values is the diameter)
ds_out = xr.Dataset({parameter: (['y', 'x'],Val.Band1.values),
                    parameter+'_sd': (['y', 'x'],Val_sd.Band1.values),
                    'lat': (['y', 'x'],coo_lat),
                    'lon': (['y', 'x'],coo_lon)},
                           coords={'y': (['y'], vecy),
                                    'x': (['x'], vecx)})

ds_out=ds_out.to_dataframe().reset_index()  #Dataframe to add the fraction of species

#Read the observed fraction of PFTs 4,5,6---------------------------------------------------------------------------

#Associate the species to a PFT as defined in the ORCHIDEE LSM. For instance, the species 1,2,3 are included in PFT 4
#Species 8 is undefinite
group_pft={"4":[1,2,3],"5":[4],"6":[5,6,7],"16":[8]}
for isp in range(1,9):
 path="ForestInventory/Tree_Species_Group/tsg_"+str(isp)+"_perc.tif"
 os.system("rm -f Species.nc")
 gdal.Translate("Species_"+str(isp)+".nc",path,format="NetCDF")
 Species=xr.open_dataset("Species_"+str(isp)+".nc")

 #Find the longitude and latitude
 x1,y1 = Species.x.values,Species.y.values

 coo_lat=np.zeros((len(y1),len(x1)))
 coo_lon=np.zeros((len(y1),len(x1)))

 for ix,xx in enumerate(x1):
  for iy,yy in enumerate(y1):

   latlong = transform.TransformPoint(xx,yy)
   coo_lat[iy,ix]=latlong[0]
   coo_lon[iy,ix]=latlong[1]


 vecx=[]; vecy=[]
 for xx in x1:
  yy=y1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecx.append(latlong[1])

 for yy in y1:
  xx=x1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecy.append(latlong[0])

 vecy=np.round(np.asarray(vecy),2)
 vecx=np.round(np.asarray(vecx),2)


 tmp = xr.Dataset({'Sp_'+str(isp): (['y', 'x'],Species.Band1.values),
                    'lat': (['y', 'x'],coo_lat),
                    'lon': (['y', 'x'],coo_lon)},
                           coords={'y': (['y'], vecy),
                                    'x': (['x'], vecx)})
 tmp=tmp.to_dataframe().reset_index()
 ds_out['Sp_'+str(isp)]=tmp['Sp_'+str(isp)].copy(deep=True)

#Associate the species to a PFT as defined in the ORCHIDEE LSM. For instance, the species 1,2,3 are included in PFT 4
for PFT in group_pft:
 ds_out["PFT_"+PFT]=0
 for ii in group_pft[PFT]:
  ds_out["PFT_"+PFT]+=ds_out["Sp_"+str(ii)] 
  del ds_out["Sp_"+str(ii)]
ds_out.to_xarray().to_netcdf(parameter.upper()+".nc")


#Read the fraction of age class -----------------------------------------------------------------------------------------

#Age classes are defined as following:
#1 0-20 years
#2 21-40 years
#3 41-60 years
#4 61-80 years
#5 81-100 years
#6 101-120 years
#7 121-140 years
#8 >140 years
#The file agecl.tif contains the dominant/most frequent age class so it should only contain values between 1 and 8.

parameter="Age"
for isp in range(1,9):
 path="ForestInventory/"+parameter+"/agecl_"+str(isp)+"_perc.tif"
 os.system("rm -f "+parameter+"_"+str(isp)+".nc")
 gdal.Translate(parameter+"_"+str(isp)+".nc",path,format="NetCDF")
 Species=xr.open_dataset(parameter+"_"+str(isp)+".nc")

 #Find the longitude and latitude
 x1,y1 = Species.x.values,Species.y.values

 coo_lat=np.zeros((len(y1),len(x1)))
 coo_lon=np.zeros((len(y1),len(x1)))

 for ix,xx in enumerate(x1):
  for iy,yy in enumerate(y1):

   latlong = transform.TransformPoint(xx,yy)
   coo_lat[iy,ix]=latlong[0]
   coo_lon[iy,ix]=latlong[1]


 vecx=[]; vecy=[]
 for xx in x1:
  yy=y1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecx.append(latlong[1])

 for yy in y1:
  xx=x1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecy.append(latlong[0])

 vecy=np.round(np.asarray(vecy),2)
 vecx=np.round(np.asarray(vecx),2)


 tmp = xr.Dataset({parameter+'_'+str(isp): (['y', 'x'],Species.Band1.values),
                    'lat': (['y', 'x'],coo_lat),
                    'lon': (['y', 'x'],coo_lon)},
                           coords={'y': (['y'], vecy),
                                    'x': (['x'], vecx)})
 tmp=tmp.to_dataframe().reset_index()
 ds_out[parameter+'_'+str(isp)]=tmp[parameter+'_'+str(isp)].copy(deep=True)


#ds_out.to_xarray().to_netcdf(parameter+'_'+str(isp)+".nc")

#Read the most dominant age class-----------------------------------------------------------------------------------
#The file agecl.tif contains the dominant/most frequent age class so it should only contain values between 1 and 8.
path="ForestInventory/"+parameter+"/agecl.tif"
os.system("rm -f "+parameter+".nc")
gdal.Translate(parameter+".nc",path,format="NetCDF")
Species=xr.open_dataset(parameter+".nc")

#Find the longitude and latitude
x1,y1 = Species.x.values,Species.y.values
coo_lat=np.zeros((len(y1),len(x1)))
coo_lon=np.zeros((len(y1),len(x1)))

for ix,xx in enumerate(x1):
  for iy,yy in enumerate(y1):

   latlong = transform.TransformPoint(xx,yy)
   coo_lat[iy,ix]=latlong[0]
   coo_lon[iy,ix]=latlong[1]


vecx=[]; vecy=[]
for xx in x1:
  yy=y1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecx.append(latlong[1])

for yy in y1:
  xx=x1[0]
  latlong = transform.TransformPoint(xx,yy)
  vecy.append(latlong[0])

vecy=np.round(np.asarray(vecy),2)
vecx=np.round(np.asarray(vecx),2)


tmp = xr.Dataset({parameter: (['y', 'x'],Species.Band1.values),
                    'lat': (['y', 'x'],coo_lat),
                    'lon': (['y', 'x'],coo_lon)},
                           coords={'y': (['y'], vecy),
                                    'x': (['x'], vecx)})
tmp=tmp.to_dataframe().reset_index()
ds_out[parameter]=tmp[parameter].copy(deep=True)
#Remove nan
ds_out=ds_out.dropna(subset=["Diameter"])
ds_out.to_xarray().to_netcdf("Diameter.nc")




