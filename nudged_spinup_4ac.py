"""

Author: Marine Remaud
Date  : March 2022
The script produces output and restart files simulating that have realistic tree diameters over European forests. The script is specific to ORCHIDEE configurated to simulate 4 diameter classes. In this case, the script optimizes the fraction of diameter class for each PFTs. The script executes the steps below:
  1) Performing a nearest neighbour interpolation of the observions from the Pucher Inventory. The interpolation is from the original 8x8 km grid to the ORCHIDEE grid. 
  2) Computation of the fraction of diameter classes within each of the grid cell
  3) Loop over the simulated years, latitudes, longitures and PFTs as defined in the ORCHIDEE LSM. For each year, grid point and PFT, the script computes the squared difference between the simulated and observed variable. The numbers are saved in a matrix J(number of years, npft,nlat,nlon)
  4) Selection of the ORCHIDEE year with a mean simulated diameter that is the closest to the observed diameter over the grid cell. The selected years are saved in a matrix J_year(npft,nlat,nlon).
  5) Loop over the latitudes, longitures and PFTs as defined in the ORCHIDEE LSM. For each ORCHIDEE grid cell and PFT, the script picks up all variables of the selected year to create a new restart file nudged toward the observations from the Pucher inventory.  

Inputs:  
 * ORCHIDEE simulations with a balanced carbon budget (very weak NBP) after a forest clear-cut
 * A file DIAMETER.nc: Observed diameter at each observed location (Pucher inventory)

Outputs: 
 * sechiba_4ac_f.nc: nudged restart file for sechiba with several  diameter classes
 * stomate_4ac_f.nc: nudged restart file for stomate with several diameter classes
 * stomate_4ac_i.nc: restart file at the end of the spinup simulation (without nudging) with several  diameter classes
 * output_4ac_f.nc : nudged ORCHIDEE output file with several diameter classes
 * output_4ac_i.nc ORCHIDEE output file at the end of the spinup simulation (without nudging) with several diameter classes

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

#---PARAMETERS------------------------------------------------------------------------------------------
list_PFT=[2,3,4,5,6,7,8,9]                                      #List of the forest PFTs
homedir="/ccc/work/cont003/gen6328/p24remau/PYTHON/SPINUP/"     #Home directory
begy=2162                                                       #Beginning of the simulation
endy=2314                                                       #End of the simulation
exp='spinup_4ac'                                                #Name of the simulation
exp2="preanspin"                                                #NAme of the job (simulation)
#Directory of the ORCHIDEE output 
dirout="/ccc/scratch/cont003/gen2201/p24remau/IGCM_OUT/OL2/TEST/"+exp+"/"+exp2+"/SBG/Output/YE/"
nb_ac=4                                                         # Number of diameter classes in ORCHIDEE
list_var=["DIAMETER"]                                        #Name of the variable to nudge
#-------------------------------------------------------------------------------------------------------

#Get the ORCHIDEE coordinates
period=str(begy)+"0101_"+str(begy)+"1231"
nf=xr.open_dataset(dirout+"/"+exp2+"_"+period+"_1Y_stomate_history.nc")
lon_orc=nf.lon.values
lat_orc=nf.lat.values
nlat=nf.dims['lat']; nlon=nf.dims['lon']; npft=nf.dims['veget']

#---NEAREST NEIGHBOUR INTERPOLATION FROM THE 8x8 km RESOLUTION TO THE ORCHIDEE RESOLUTION----------------
Obs=xr.open_dataset(list_var[0]+".nc").to_dataframe()                          #Read observed diameter
Obs.rename(columns={"Diameter":"obs","Diameter_sd":"error"},inplace=True)
Obs["Parameter"]=list_var[0]
Obs["clon"]=Obs.apply(lambda row: np.abs(lon_orc-row.lon).argmin(),axis=1) #Find the nearest longitude for each observation
Obs["clat"]=Obs.apply(lambda row: np.abs(lat_orc-row.lat).argmin(),axis=1) #Find the nearest latiude (in ORCHIDEE) for each observation
Obs1=Obs.groupby(["clon","clat"]).mean().reset_index()                     #Average diameters over ORCHIDEE grid cell (used as a template)

#---COMPUTATION OF THE FRACTION OF DIAMETER CLASSES PER ORCHIDEE GRID CELL------------------------------- 
#Declaration of columns C1,..,C4 :fraction of age class as defined in the ORCHIDEE LSM
Obs1["C1"]=0; Obs1["C2"]=0; Obs1["C3"]=0; Obs1["C4"]=0

#Computation of fraction of age classes
for irr,row in Obs1.iterrows():
 mask=(Obs.clon==int(row.clon))&(Obs.clat==int(row.clat))&(Obs.Parameter=="DIAMETER").copy(deep=True)
 Obs1.loc[irr,"C1"]=float(len(Obs[mask&(Obs.obs<7)]))               /float(len(Obs[mask]))
 Obs1.loc[irr,"C2"]=float(len(Obs[mask&(Obs.obs>=7)&(Obs.obs<20)])) /float(len(Obs[mask]))
 Obs1.loc[irr,"C3"]=float(len(Obs[mask&(Obs.obs>=20)&(Obs.obs<40)]))/float(len(Obs[mask]))
 Obs1.loc[irr,"C4"]=float(len(Obs[mask&(Obs.obs>=40)]))             /float(len(Obs[mask]))
Obs1.dropna(inplace=True)

#---Computation of the square of the difference between the observed and simulated values----------------
J=np.zeros((nlat,nlon,npft,(endy-begy+1)))                       #J(j,i,k,year)=(diameter_obs-diameter_sim(year))^2
J[:,:,:,:]=np.nan
J_year=np.zeros((nlat,nlon,npft))                                #Matrix filled with year for which the differences between the simulated ans observed values are minimal: J_year(j,i,k)=min(J(j,i,k,:))                               
for iv,var in enumerate(list_var):
 for yy in range(begy,endy+1):                                   #Loop over the years of ORCHIDEE spinup simulattion
  period=str(yy)+"0101_"+str(yy)+"1231"                          #Opening of the ORCHIDEE yearly output for the year yy
  namef=exp2+"_"+period+"_1Y_stomate_history.nc"
  f = Dataset(dirout+"/"+namef)
  ncf=f.variables["VEGET_MAX"][-1,:,:,:]                         #The fraction of PFT is fixed and taken from the ORCHIDEE configuration (ESA maps are prescribed)
  f.close()
  #Loop over ORCHIDEE grid cells, PFTs and diameter classes
  for ilat in range(nlat):
   for ilon in range(nlon):
    for ipft in list_PFT:                                         
     for ac in range(1,nb_ac+1):               
      num_PFT=1+(ipft-2)*4+ac
      Obs1_tmp=Obs1[(Obs1.clon==ilon)&(Obs1.clat==ilat)].copy(deep=True)
      if Obs1_tmp.empty: continue
      #Extract the grid cell closest to the observation location 
      ncf_ij=np.squeeze(ncf.data[num_PFT,ilat,ilon])
      i_beg=1+(ipft-2)*4
      i_end=1+(ipft-2)*4+nb_ac+1
      veget_frac=np.sum(ncf[i_beg:i_end,ilat,ilon].data)         #Fraction of the PFT within  the grid cell (ilat,ilon)
      #Square of the difference between the observed and simulated values
      error=np.abs(ncf_ij-Obs1_tmp["C"+str(ac)].iloc[0]*veget_frac)
      if np.isnan(J[ilat,ilon,num_PFT-1,yy-begy]): 
       J[ilat,ilon,num_PFT-1,yy-begy]=0
      J[ilat,ilon,num_PFT-1,yy-begy] =error

#Select the year for each pixel and each PFT
id_nan=np.where(~np.isnan(J))
id_lat=np.unique(id_nan[0]); id_lon=np.unique(id_nan[1]);id_pft=np.unique(id_nan[2])
for ilat in id_lat:
 for ilon in id_lon:
  for ipft in id_pft:
   if len(J[ilat,ilon,ipft,~np.isnan(J[ilat,ilon,ipft,:])] )==0: continue
   J_year[ilat,ilon,ipft]=begy+np.nanargmin(J[ilat,ilon,ipft,:], axis=-1)

###Year to choose when nan
year_fill_J=Counter( J_year[J_year!=0].flatten()).most_common(1)[0][0]
J_year[J_year==0]=year_fill_J

#CREATION OF THE RESTART FILE------------------------------------------------------------------------------------------
name_restart={"SRF":"sechiba","SBG":"stomate"}
for rr in name_restart.keys(): 
 #Load a restart file as a template
 print(dirout+"/../../../"+rr+"/Restart/"+exp2+"_"+str(begy)+"1231_"+name_restart[rr]+"_rest.nc")
 file_restart_0=dirout+"/../../../"+rr+"/Restart/"+exp2+"_"+str(begy)+"1231_"+name_restart[rr]+"_rest.nc"
 restart=xr.open_dataset(file_restart_0,decode_times=False,decode_cf=False)
 for yy in range(begy,endy+1):
  file_restart=dirout+"/../../../"+rr+"/Restart/"+exp2+"_"+str(yy)+"1231_"+name_restart[rr]+"_rest.nc"
  id_yr=np.where(J_year==yy) 
  id_lon=id_yr[1]
  id_lat=id_yr[0]
  id_pft=id_yr[2]
  if len(id_yr)==0: continue
  for var in restart.keys():
   f = Dataset(file_restart)
   ncf=f.variables[var][:]
   f.close()
   dim_var=np.shape(ncf)
   if (len(dim_var)==1)|(len(ncf)==1): continue
   loc_lat=np.where(np.asarray(dim_var)==nlon)[0][0]
   loc_pft=np.where(np.asarray(dim_var)==npft)[0]
   print(loc_pft,len(dim_var),var)
   if loc_pft == 0: continue
#   print(dim_var,var,yy)   
   if (loc_pft==1)&(len(dim_var)==4):
    restart[var].values[:,id_pft,id_lat,id_lon]=ncf[:,id_pft,id_lat,id_lon]
    print("selection")
   elif (loc_pft==2)&(len(dim_var)==5):
    restart[var].values[:,:,id_pft,id_lat,id_lon]=ncf[:,:,id_pft,id_lat,id_lon]
   elif (loc_pft==3)&(len(dim_var)==6):
    restart[var].values[:,:,:,id_pft,id_lat,id_lon]=ncf[:,:,:,id_pft,id_lat,id_lon]
 os.system("rm -f "+homedir+"/"+name_restart[rr]+"_4ac_f.nc")
 os.system("rm -f "+homedir+"/"+name_restart[rr]+"_4ac_i.nc")
 
 restart.to_netcdf(homedir+"/"+name_restart[rr]+"_4ac_f.nc") 
 os.system("cp -r "+file_restart+ " "+homedir+"/"+name_restart[rr]+"_4ac_i.nc") 
 

#CREATION OF THE OUTPUT FILE------------------------------------------------------------------------------------------
print('Creation of the output files')
#Load an ORCHIDEE outputt file as a template
file_output_0=dirout+exp2+"_"+period+"_1Y_stomate_history.nc"
restart=xr.open_dataset(file_output_0,decode_times=False,decode_cf=False)
list_vars=["HEIGHT","DIAMETER","VEGET_MAX"]
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
   elif (loc_pft==2)&(len(dim_var)==5):
    restart[var].values[:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,id_pft,id_lat,id_lon])
   elif (loc_pft==3)&(len(dim_var)==6):
    restart[var].values[:,:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,:,id_pft,id_lat,id_lon])
os.system("rm -f "+homedir+"/output_4ac_f.nc")
os.system("rm -f "+homedir+"/output_4ac_i.nc")
    
restart.to_netcdf(homedir+"/output_4ac_f.nc")
os.system("cp -r "+file_output+ " "+homedir+"/output_4ac_i.nc")


