"""

Author: Marine Remaud
Date  : March 2022
The script produces output and restart files that have realistic tree diameters over European forests. The script executes the steps below:
  1) Perform a nearest neighbour interpolation of the observed diameter from the Pucher Inventory. The interpolation is from the original 8x8 km grid to the ORCHIDEE grid. 
  2) Loop over the simulated years, latitudes, longitures and PFTs as defined in the ORCHIDEE LSM. For each year, grid point and PFT, the script computes the squared difference between the simulated and observed variable. The numbers are saved in a matrix J(number of years, npft,nlat,nlon)
  3) Selection of the ORCHIDEE year with a mean simulated diameter that is the closest to the observed diameter over the grid cell. The selected years are saved in a matrix J_year(npft,nlat,nlon).
  4) Loop over the latitudes, longitures and PFTs as defined in the ORCHIDEE LSM. For each ORCHIDEE grid cell and PFT, the script picks up all variables of the selected year to create a new restart file nudged toward the observations from the Pucher inventory.  

Inputs:  
 * ORCHIDEE simulations with a balanced carbon budget (very weak NBP) after a forest clear-cut
 * A file DIAMETER.nc with the observed diameter from the Pucher inventory

Outputs: 
 * sechiba_f.nc: nudged restart file for sechiba 
 * stomate_f.nc: nudged restart file for stomate
 * stomate_i.nc: restart file at the end of the spinup simulation (without nudging)
 * output_f.nc : nudged ORCHIDEE output file 
 * output_i.nc ORCHIDEE output file at the end of the spinup simulation (without nudging)
"""

import pandas as pd
import numpy as np
import os
import xarray as xr
from netCDF4 import Dataset
import datetime
import matplotlib.pyplot as plt
import copy
import calendar
from collections import Counter


#LIST OF PARAMETERS---------------------------------------------------------------------------------------------------------
homedir="/ccc/work/cont003/gen6328/p24remau/PYTHON/SPINUP/" #Homedir
begy=2162                                                   #first year of the ORCHIDEE spinup simulation after the clearcut
endy=2314                                                   #last year of the ORCHIDEE spinup simulation
nb_PFT=15                                                   # Number of PFT as defined in the spin-up ORCHIDEE configuration
exp="spinup3_n" #'''spinup_3006'                            # Name 1 of the simulation  
exp2="anspin"                                               # Name 2 of the simulation
dirout="/ccc/scratch/cont003/gen2201/p24remau/IGCM_OUT/OL2/TEST/"+exp+"/"+exp2+"/SBG/Output/YE/"   #Path toward the ORCHIDEE output
list_var=["DIAMETER"]                                       #Name of the observed variable 
list_pft=[4,5,6,7,8,9]                                      #PFTs covered by forest

#READ OBSERVATIONS---------------------------------------------------------------------------------------------------------
Obs1=pd.DataFrame()
for iv,var in enumerate(list_var):
 tmp=xr.open_dataset(var+".nc").to_dataframe()
 tmp.rename(columns={"Diameter":"obs","Diameter_sd":"error"},inplace=True)
 tmp["Parameter"]=var
 Obs1=Obs1.append(tmp)
 Obs1["obs"]*=10**(-2)  #Conversion from cm to m
 Obs1["error"]*=10**(-2)

#Get the ORCHIDEE coordinates----------------------------------------------------------------------------------------------
period=str(begy)+"0101_"+str(begy)+"1231"
nf=xr.open_dataset(dirout+"/"+exp2+"_"+period+"_1Y_stomate_history.nc")
lon_orc=nf.lon.values
lat_orc=nf.lat.values
nlat=nf.dims['lat']; nlon=nf.dims['lon']; npft=nf.dims['veget']


#NEAREST NEIGHBOUR INTERPOLATION---------------------------------------------------------------------------------------------
Obs1["clon"]=Obs1.apply(lambda row: np.abs(lon_orc-row.lon).argmin(),axis=1) #Find the closest ORCHIDEE longitude clon to each observation
Obs1["clat"]=Obs1.apply(lambda row: np.abs(lat_orc-row.lat).argmin(),axis=1) #Find the closest ORCHIDEE latitude clat to each observation
Obs1.dropna(inplace=True)
#Obs1=Obs1.groupby(["clon""clat"]).mean().reset_index() #Average over each ORCHIDEE grid cell

#Boreal PFTs
#Obs1["PFT_7"]=0; Obs1["PFT_8"]=0; Obs1["PFT_9"]=0
#Obs1["PFT_7"]=Obs1.apply(lambda row: row.PFT_4 if nf.VEGET_MAX.values[0,3,row.clat,row.clon]>0 else 0,axis=1)
#Obs1["PFT_8"]=Obs1.apply(lambda row: row.PFT_5 if nf.VEGET_MAX.values[0,4,row.clat,row.clon]>0 else 0,axis=1)
#Obs1["PFT_9"]=Obs1.apply(lambda row: row.PFT_6 if nf.VEGET_MAX.values[0,5,row.clat,row.clon]>0 else 0,axis=1)
#Obs1["PFT_4"]=Obs1.apply(lambda row: 0 if  row.PFT_7>0 else row.PFT_4,axis=1)
#Obs1["PFT_5"]=Obs1.apply(lambda row: 0 if  row.PFT_8>0 else row.PFT_5,axis=1)
#Obs1["PFT_6"]=Obs1.apply(lambda row: 0 if  row.PFT_9>0 else row.PFT_6 ,axis=1)
#Obs1["PFT"]=Obs1.apply(lambda row: 4+np.argmax([row.PFT_4,row.PFT_5,row.PFT_6,row.PFT_7,row.PFT_8,row.PFT_9]),axis=1)
#Obs1["PFT"]=Obs1.apply(lambda row: int(row.PFT) if row["PFT_"+str(row.PFT)]>0.5 else np.nan ,axis=1)
#Obs1["PFT"]=Obs1.apply(lambda row: int(row.PFT),axis=1)

#Computation of the square of the difference between the observed and simulated values------------------------------------------

J=np.zeros((nlat,nlon,npft,(endy-begy+1)))      #J(j,i,k,year)=(diameter_obs-diameter_sim(year))^2
J[:,:,:,:]=np.nan                               
Jmin=np.zeros((nlat,nlon,npft,(endy-begy+1)))   #Jmin(j,i,k,year)=(diameter_obs-error_obs-diameter_sim(year))^2
Jmin[:,:,:,:]=np.nan
Jmax=np.zeros((nlat,nlon,npft,(endy-begy+1)))   #Jmax(j,i,k,year)(diameter_obs+error_obs-diameter_sim(year))^2
Jmax[:,:,:,:]=np.nan

J_year=np.zeros((nlat,nlon,npft))               #Matrix filled with year for which the differences between the simulated ans observed values are minimal: J_year(j,i,k)=min(J(j,i,k,:))
Jmin_year=np.zeros((nlat,nlon,npft))            #Jmin_year(j,i,k)=min(Jmin(j,i,k,:)
Jmax_year=np.zeros((nlat,nlon,npft))            #Jmax_year(j,i,k)=min(Jmax(j,i,k,:)


debug=np.zeros((1,6,endy-begy+1))
for iv,var in enumerate(list_var):
 for yy in range(begy,endy+1):
  Obs1_yr=Obs1[(Obs1.Parameter==var)].copy(deep=True)
  if Obs1_yr.empty: continue
  #Open the ORCHIDEE output file for the year yy to obtain the simulated parameter (e.g diameter)
  period=str(yy)+"0101_"+str(yy)+"1231"
  namef=exp2+"_"+period+"_1Y_stomate_history.nc"
  f = Dataset(dirout+"/"+namef)
  ncf=f.variables[var][-1,:,:,:]
  f.close()
  for ilat in Obs1_yr.clat.unique():     #Loop over the ORCHIDEE latitudes
   for ilon in Obs1_yr[Obs1_yr.clat==ilat].clon.unique():  #Loop over the ORCHIDEE longitudes
    for ipft in list_pft:
     if (ilat==13)&(ilon==11): debug[0,:,yy-begy]=np.squeeze(ncf[3:9,ilat,ilon])   #Debug
     Obs1_tmp=Obs1_yr[(Obs1_yr.clon==ilon)&(Obs1_yr.clat==ilat)].copy(deep=True)
     if Obs1.empty: continue
     #Extract the grid cell at the observed location 
     ncf_ij=np.squeeze(ncf[int(ipft-1),ilat,ilon])
     if np.isnan(J[ilat,ilon,ipft-1,yy-begy]): 
      J[ilat,ilon,ipft-1,yy-begy]=0
      Jmin[ilat,ilon,ipft-1,yy-begy]=0
      Jmax[ilat,ilon,ipft-1,yy-begy]=0
     #Square of the difference between the observed and simulated value
     J[ilat,ilon,ipft-1,yy-begy]=(ncf_ij-Obs1_tmp["obs"].mean())**2
     Jmin[ilat,ilon,ipft-1,yy-begy]=(ncf_ij-Obs1_tmp["obs"].mean()-Obs1_tmp["error"].std())**2
     Jmax[ilat,ilon,ipft-1,yy-begy]=(ncf_ij-Obs1_tmp["obs"].mean()+Obs1_tmp["error"].std())**2

#Selection of the simulated year in ORCHIDEE for each pixel and PFT---------------------------------------------------
#The year minimizes the diffence between the simulated and observed value of the chosen parameter (e.g. diameter)
id_nan=np.where(~np.isnan(J))
id_lat=np.unique(id_nan[0]); id_lon=np.unique(id_nan[1]);id_pft=np.unique(id_nan[2])
for ilat in id_lat:
 for ilon in id_lon:
  for ipft in id_pft:
   if len(J[ilat,ilon,ipft,~np.isnan(J[ilat,ilon,ipft,:])] )==0: continue
   J_year[ilat,ilon,ipft]=begy+np.nanargmin(J[ilat,ilon,ipft,:], axis=-1)
   Jmin_year[ilat,ilon,ipft]=begy+np.nanargmin(Jmin[ilat,ilon,ipft,:], axis=-1)
   Jmax_year[ilat,ilon,ipft]=begy+np.nanargmin(Jmax[ilat,ilon,ipft,:], axis=-1)

#Year to choose when nan
year_fill_J=Counter( J_year[J_year!=0].flatten()).most_common(1)[0][0]
year_fill_Jmin=Counter( Jmin_year[Jmin_year!=0].flatten()).most_common(1)[0][0]
year_fill_Jmax=Counter( Jmax_year[Jmax_year!=0].flatten()).most_common(1)[0][0]
J_year[J_year==0]=year_fill_J
Jmax_year[Jmax_year==0]=year_fill_Jmax
Jmin_year[Jmin_year==0]=year_fill_Jmin

#CREATION OF THE RESTART FILE------------------------------------------------------------------------------------------

print('Creation of the restart file')
name_restart={"SRF":"sechiba","SBG":"stomate"}

for rr in name_restart.keys(): #loop over sechiba, stomate
 #Load a restart file as a template
 file_restart_0=dirout+"/../../../"+rr+"/Restart/"+exp2+"_"+str(begy)+"1231_"+name_restart[rr]+"_rest.nc"
 restart=xr.open_dataset(file_restart_0,decode_times=False,decode_cf=False)
 for yy in range(begy,endy+1): #Loop over the spinup years
  file_restart=dirout+"/../../../"+rr+"/Restart/"+exp2+"_"+str(yy)+"1231_"+name_restart[rr]+"_rest.nc"
  id_yr=np.where(J_year==yy) 
  id_lon=id_yr[1]
  id_lat=id_yr[0]
  id_pft=id_yr[2]
  if len(id_yr)==0: continue
  for var in restart.keys():
   #Open the restart file for the spinup year yy 
   f = Dataset(file_restart)
   ncf=f.variables[var][:]
   f.close()
   dim_var=np.shape(ncf)
   if (len(dim_var)==1)|(len(ncf)==1): continue
   loc_lat=np.where(np.asarray(dim_var)==nlon)[0][0]
   loc_pft=np.where(np.asarray(dim_var)==npft)[0]
   if loc_pft == 0: continue
   if (loc_pft==1)&(len(dim_var)==4):
    restart[var].values[:,id_pft,id_lat,id_lon]=np.copy(ncf[:,id_pft,id_lat,id_lon])
   elif (loc_pft==2)&(len(dim_var)==5):
    restart[var].values[:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,id_pft,id_lat,id_lon])
   elif (loc_pft==3)&(len(dim_var)==6):
    restart[var].values[:,:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,:,id_pft,id_lat,id_lon])
 os.system("rm -f "+homedir+"/"+name_restart[rr]+"_1ac_f.nc")
 os.system("rm -f "+homedir+"/"+name_restart[rr]+"_1ac_i.nc")
 
 restart.to_netcdf(homedir+"/"+name_restart[rr]+"_1ac_f.nc") 
 os.system("cp -r "+file_restart+ " "+homedir+"/"+name_restart[rr]+"_1ac_i.nc") 

#CREATION OF THE OUTPUT FILE------------------------------------------------------------------------------------------
print('Creation of the output files')
#Load an output file as a template#####
file_output_0=dirout+exp2+"_"+period+"_1Y_stomate_history.nc"
restart=xr.open_dataset(file_output_0,decode_times=False,decode_cf=False)
list_vars=["HEIGHT","DIAMETER"]
for yy in range(begy,endy+1):
  period=str(yy)+"0101_"+str(yy)+"1231"
  file_output=dirout+exp2+"_"+period+"_1Y_stomate_history.nc"
  id_yr=np.where(J_year==yy)
  id_lon=id_yr[1]
  id_lat=id_yr[0]
  id_pft=id_yr[2]
  if len(id_yr)==0: continue
  for var in list_vars:
   f = Dataset(file_output)
   ncf=f.variables[var][:]
   f.close()
   dim_var=np.shape(ncf)
   if (len(dim_var)==1)|(ncf.count()==1): continue
   loc_lat=np.where(np.asarray(dim_var)==nlon)[0][0]
   loc_pft=np.where(np.asarray(dim_var)==npft)[0][0]
   if loc_pft == 0: continue
   if (loc_pft==1)&(len(dim_var)==4):
    restart[var].values[:,id_pft,id_lat,id_lon]=np.copy(ncf[:,id_pft,id_lat,id_lon])
    print("selection")
   elif (loc_pft==2)&(len(dim_var)==5):
    restart[var].values[:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,id_pft,id_lat,id_lon])
   elif (loc_pft==3)&(len(dim_var)==6):
    restart[var].values[:,:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,:,id_pft,id_lat,id_lon])
    print(var)
os.system("rm -f "+homedir+"/output_1ac_f.nc")
os.system("rm -f "+homedir+"/output_1ac_i.nc")
    
restart.to_netcdf(homedir+"/output_1ac_f.nc")
os.system("cp -r "+file_output+ " "+homedir+"/output_1ac_i.nc")
 
