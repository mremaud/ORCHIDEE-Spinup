"""

Author: Marine Remaud
Date: Mai 2022
European maps of the variables obtained before and after the nudging. 
For instance for the diameter, the script compares the diameter at the end of the spinup and the diameter that has been nudged toward the observations.

Requirements: Perform a spinup simulation with ORCHIDEE beforehand. The files stomate_i.nc, output_i.nc are the restart and the output file respectively at the end of the spiinup. 
Inputs:
  - stomate_nac_i.nc: Restart file at the end of the spinup with n diameter classes
  - output_nac_i.nc: Output file (SRF) at the end of the spinup with n diameter classes
  - stomate_nac_f.nc: Nudged restart file with n diameter classes
  - output_nac_f.nc: Nudged output file with n diameter classes

"""

import matplotlib as mpl
import xarray as xr
import pandas as pd
import numpy as np
import os
from netCDF4 import Dataset
import datetime
import copy
import calendar
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as mcolors

#---PARAMETERS------------------------------------------------------------------------------------------------------------------------
begy=2162           #First year of the spinup ORCHIDEE simulation after the forest clearcut
endy=2314           #End year
exp='spinup_3006'   #Name of the simulation
exp2="preanspin"    #Name of the job
dirout="/ccc/scratch/cont003/gen2201/p24remau/IGCM_OUT/OL2/TEST/"+exp+"/"+exp2+"/SBG/Output/YE/"  #Path toward the ORCHIDEE simulation
nb_ac=4             #Number of age classes
#Coordonate of the map 
latmin=35
latmax=60
lonmin=-9
lonmax=29
list_var_restart=["soil_carbon","litter_c"] #,"resp_hetero"]   #List of the variables to plot in the restart file
#--------------------------------------------------------------------------------------------------------------------------------------

#Get the longitude, latitude coordinates
period=str(begy)+"0101_"+str(begy)+"1231"
nf=xr.open_dataset(dirout+"/"+exp2+"_"+period+"_1Y_stomate_history.nc")
lon_orc=nf.lon.values
lat_orc=nf.lat.values
nlat=nf.dims['lat']; nlon=nf.dims['lon']; npft=nf.dims['veget']

#Relative difference between the variables in the restart file (e.g. soil_carbon) -------------------------------- 
for var in list_var_restart: 
 prior=xr.open_dataset("stomate_"+str(nb_ac)+"ac_i.nc",decode_times=False)[var].values
 post=xr.open_dataset("stomate_"+str(nb_ac)+"ac_f.nc",decode_times=False)[var].values
 prior=np.squeeze(prior)
 post=np.squeeze(post)
 print(np.shape(post))
 if len(np.shape(prior)) > 2: 
  prior=np.squeeze(np.mean(prior,axis=0))
 if len(np.shape(prior))==4: prior=np.mean(prior,axis=0)
 if len(np.shape(prior))==3: prior=np.mean(prior,axis=0)
 prior[prior== 1.00000000e+20]=np.nan

 if len(np.shape(post)) > 2:
  post=np.squeeze(np.mean(post,axis=0))
 if len(np.shape(post))==4: post=np.squeeze(np.mean(post,axis=0))
 if len(np.shape(post))==3: post=np.mean(post,axis=0)
 post[post== 1.00000000e+20]=np.nan

 diff=post-prior
 #Map showing the impact of the nudging
 fig, ax = plt.subplots()
 X,Y=np.meshgrid(lon_orc,lat_orc)
 m = Basemap(projection='cyl',resolution='c',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
 m.drawcoastlines()
 m.contourf(X,Y,diff/prior,cmap="RdBu_r",levels=np.linspace(-1.1, 1.1, 12))
 #m.contourf(X,Y,prior,cmap="Spectral_r")

 plt.colorbar()
 plt.title(var)
 n_graticules = 18
 parallels = np.arange(-90., 90, 20)
 meridians = np.arange(0., 360., 40)
 lw = 0.5
 dashes = [1,1] # 5 dots, 7 spaces... repeat
 graticules_color = 'grey'

#----------------------------------------------------------------------------------------------------------------
#Evaluation of the nudging processes: comparison between the nudged and the observed diameter 

list_var=["HEIGHT","DIAMETER"]
list_var=["DIAMETER"]
PFTi=2+nb_ac*2
PFTe=1+nb_ac*8+1
list_pft=np.arange(PFTi,PFTe)
print(list_pft)
for var in list_var:
 prior=xr.open_dataset("output_"+str(nb_ac)+"ac_i.nc",decode_times=False,decode_cf=False)
 post=xr.open_dataset("output_"+str(nb_ac)+"ac_f.nc",decode_times=False,decode_cf=False)
 prior=prior.sel(veget=list_pft)   #Select the forested PFTs
 post=post.sel(veget=list_pft)     
 if (var == "DIAMETER")|(var=="HEIGHT"):
  prior=(prior[var]*prior["VEGET_MAX"]/prior["VEGET_MAX"].sum("veget"))  #Weighted over the vegetation fraction
  post=(post[var]*post["VEGET_MAX"]/post["VEGET_MAX"].sum("veget"))
  prior=prior.sum("veget").values
  post =post.sum("veget").values

 prior=np.squeeze(prior)
 post=np.squeeze(post)
 if len(np.shape(prior))==3: prior=np.mean(prior,axis=0)
 if len(np.shape(post))==3: post=np.mean(post,axis=0)
 post[post== 1.00000000e+20]=np.nan
 prior[prior== 1.00000000e+20]=np.nan

 diff=post #-prior/prior  #Impact of the nudging
 diff=np.flip(diff,axis=0)*10**2

 fig, ax = plt.subplots()
 X,Y=np.meshgrid(lon_orc-1,np.flip(lat_orc,axis=0)-1)
 m = Basemap(projection='cyl',resolution='c',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
 m.drawcoastlines()
 post[post>10**10]=np.nan
 prior[prior>10**10]=np.nan
 m.pcolormesh(X,Y,diff,vmin=9,vmax=32,cmap="Spectral_r") 
 plt.colorbar()
 plt.title(var)
 n_graticules = 18
 parallels = np.arange(-90., 90, 2)
 meridians = np.arange(0., 360., 2)
 lw = 0.5
 dashes = [1,1] # 5 dots, 7 spaces... repeat
 graticules_color = 'grey'
 m.drawparallels(parallels, linewidth=0.2,labels=[1,0,0,0], dashes=dashes, color=graticules_color,zorder=20)
 m.drawmeridians(meridians,linewidth=0.2, labels=[0,0,0,1],dashes=dashes, color=graticules_color,zorder=20)

#Map of the observations----------------------------------------------------------------------------------------
Diam=xr.open_dataset("DIAMETER.nc").to_dataframe()
Diam=Diam.dropna()
#Nearest neighbour interpolation
Diam["clon"]=Diam.apply(lambda row: np.abs(lon_orc-row.lon).argmin(),axis=1)
Diam["clat"]=Diam.apply(lambda row: np.abs(lat_orc-row.lat).argmin(),axis=1)
Diam["lon"]=Diam.apply(lambda row: lon_orc[int(row.clon)],axis=1)
Diam["lat"]=Diam.apply(lambda row: lat_orc[int(row.clat)],axis=1)
Diam=Diam.groupby(["clon","clat"]).mean().reset_index()

df_pv=Diam.pivot_table(index="lat", columns="lon").Diameter
da = xr.DataArray(data=df_pv)

fig, ax = plt.subplots(1,1)
m = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,\
              llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='c')
cmap=plt.cm.Spectral_r  #plt.cm.PiYG_r
m.drawcoastlines()
n_graticules = 18
parallels = np.arange(-90., 170, 2)
meridians = np.arange(0., 480., 2)
lw = 0.5; dashes = [1,1] # 5 dots, 7 spaces... repeat
graticules_color = 'grey'
m.drawparallels(parallels, linewidth=0.5,labels=[1,0,0,0], fontsize=11,dashes=dashes, color=graticules_color,zorder=20)
m.drawmeridians(meridians,linewidth=0.5, labels=[0,1,0,1],fontsize=11,dashes=dashes, color=graticules_color,zorder=20)
Y, X=np.meshgrid(da.lon-1,da.lat-1)
Y, X = m(Y,X)
ax1=m.pcolormesh(Y,X,da.values,cmap=cmap,vmin=9,vmax=32)
cb = mpl.colorbar.ColorbarBase(ax1,orientation="horizontal", spacing='proportional', extend="both", format='%1i')
cb.set_label('Diameter (cm)',fontsize=13)



