"""

Author: Marine Remaud
Date  : Mai 2023
The script produces output and restart files simulating that have realistic tree heights over northern american forests. It minimizes the difference between the observed and simulated diameters. The script executes the steps below:
  1) Performing a nearest neighbour interpolation of the observions. The interpolation is from the original grid to the ORCHIDEE grid. 
  2) Loop over the simulated years, latitudes, longitures and PFTs as defined in the ORCHIDEE LSM. For each year, grid point and PFT, the script computes the squared difference between the simulated and observed variable. The numbers are saved in a matrix J(number of years, npft,nlat,nlon)
  3) Selection of the ORCHIDEE year with a mean simulated diameter that is the closest to the observed diameter over the grid cell. The selected years are saved in a matrix J_year(npft,nlat,nlon).
  4) Loop over the latitudes, longitures and PFTs as defined in the ORCHIDEE LSM. For each ORCHIDEE grid cell and PFT, the script picks up all variables of the selected year to create a new restart file nudged toward the observed heights.  

Inputs:  
 * ORCHIDEE simulations with a balanced carbon budget (very weak NBP) after a forest clear-cut
 * A file HEIGHT_DOM.nc: Observed height at each observed location 

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

#---PARAMETERS (TO CHANGE)------------------------------------------------------------------------------------------
list_PFT=[2,3,4,5,6,7,8,9]                                              #List of the PFTs covered by forests
homedir="/home/users/mremaud/PYTHON/SPINUP/PYTHON/"                     #Home directory
dirobs="/home/users/mremaud/PYTHON/SPINUP/PYTHON/"
begy=1900                                                               #Beginning of the simulation
endy=2019                                                               #End of the simulation
exp='spinup_fm2'                                                         #Name of the simulation
exp2="anspin"               
exp2="FG2"                                            #Name of the job (simulation)
dirout="/home/scratch01/mremaud/IGCM_OUT/OL2/TEST/"+exp+"/"+exp2+"/SBG/Output/MO/"   #Directory of the ORCHIDEE output
def_file="/home/users/mremaud/ORCHIDEE/modipsl/config/ORCHIDEE_OL/TEST1D_DO/PARAM/orchidee_pft.def_15pft.1ac"   #Def file
list_var=["HEIGHT_DOM"]                                        #Name of the variable to nudge

#-------------------------------------------------------------------------------------------------------
#Read the def_file that contains the number  of age class, pft
with open(def_file) as f:
    for line in f:
     if "NAGEC" in line:
      nb_ac=int(line[6:7])
     if "NVMAP" in line:
      n_veget=int(line[6:8])
     if "NVM" in line:
      npft=int(line[6:8])

#Table_PFT= correspondance between metaclass and pft
#Table_AGE= correspondance between metaclass and age
Table_PFT=np.zeros(npft)
Table_AGE=np.zeros(npft)
with open(def_file) as f:
    for line in f:
      if "AGEC_GROUP" in line:
       i_mc=int(line[12:14])
       Table_PFT[i_mc-1]=int(line[15:17])
      if ("PFT_NAME" in line )&("age" in line):
       i_mc=int(line[10:12])
       Table_AGE[i_mc-1]=int(line[-2:])

print("Number of meta-classes",npft)
print("Number of age classes",nb_ac)
print("Number of PFTs",n_veget)
print(Table_PFT, Table_AGE)

#---------------------------------------------------------------------------------------------------------
#Get the ORCHIDEE coordinates
period=str(begy)+"0101_"+str(begy)+"1231"
nf=xr.open_dataset(dirout+"/"+exp2+"_"+period+"_1M_stomate_history.nc")
lon_orc=nf.lon.values
lat_orc=nf.lat.values
nlat=nf.dims['lat']; nlon=nf.dims['lon']; npft=nf.dims['veget']

#---NEAREST NEIGHBOUR INTERPOLATION FROM THE RESOLUTION OF THE OBSERVATIONS TO THE ORCHIDEE RESOLUTION----------------
Obs=xr.open_dataset(dirobs+list_var[0]+".nc").to_dataframe()                          #Read observed diameter
Obs["Parameter"]=list_var[0]
Obs["clon"]=Obs.apply(lambda row: np.abs(lon_orc-row.lon).argmin(),axis=1) #Find the nearest longitude for each observation
Obs["clat"]=Obs.apply(lambda row: np.abs(lat_orc-row.lat).argmin(),axis=1) #Find the nearest latiude (in ORCHIDEE) for each observation
Obs1=Obs.groupby(["clon","clat"]).mean(numeric_only=True).reset_index()                     #Average diameters over ORCHIDEE grid cell (used as a template)


#---Computation of the square of the difference between the observed and simulated values----------------
J=np.zeros((nlat,nlon,npft,(endy-begy+1)))                       #J(j,i,k,year)=(diameter_obs-diameter_sim(year))^2
J[:,:,:,:]=np.nan
J_year=np.zeros((nlat,nlon,npft))                                #Matrix filled with year for which the differences between the simulated ans observed values are minimal: J_year(j,i,k)=min(J(j,i,k,:))                               
for iv,var in enumerate(list_var):
 for yy in range(begy,endy+1):                                   #Loop over the years of ORCHIDEE spinup simulattion
  period=str(yy)+"0101_"+str(yy)+"1231"                          #Opening of the ORCHIDEE yearly output for the year yy
  namef=exp2+"_"+period+"_1M_stomate_history.nc"
  f = Dataset(dirout+"/"+namef)
  ncf_diam=f.variables[var][-1,:,:,:]                            #Diameter (4 ac)
  f.close()
  #Loop over ORCHIDEE grid cells, PFTs and diameter classes
  for ilat in range(nlat):
   for ilon in range(nlon):
    Obs1_tmp=Obs1[(Obs1.clon==ilon)&(Obs1.clat==ilat)].copy(deep=True)
    if Obs1_tmp.empty: continue
    for ipft in list_PFT:                             
     Obs1_tmp=Obs1[(Obs1.clon==ilon)&(Obs1.clat==ilat)&(Obs1.pft==ipft)].copy(deep=True)
     if Obs1_tmp.empty: continue
     ind_mc=np.where(Table_PFT==ipft)[0]
     ind_ac=Table_AGE[ind_mc]
     error = np.zeros(nb_ac)
     diam  = np.zeros(nb_ac)
     for imc in ind_mc:               
      ac=int(Table_AGE[imc])
      #Extract the grid cell closest to the observation location 
      ncf_ij=np.copy(np.squeeze(ncf_diam.data[imc,ilat,ilon]))
      diam[ac-1]=np.copy(ncf_ij)
      if np.isnan(J[ilat,ilon,imc,yy-begy]): 
       J[ilat,ilon,imc,yy-begy]=0
     if len(diam[diam!=0])!=0:  diam[diam==0]=np.max(diam)
     J[ilat,ilon,ind_mc,yy-begy]=np.min(np.abs(diam-Obs1_tmp.obs.iloc[0]))

#Select the year for each pixel and each PFT
id_nan=np.where(~np.isnan(J))
id_lat=np.unique(id_nan[0]); id_lon=np.unique(id_nan[1]);id_pft=np.unique(id_nan[2])
for ilat in id_lat:
 for ilon in id_lon:
  for ipft in id_pft:
   if len(J[ilat,ilon,ipft,~np.isnan(J[ilat,ilon,ipft,:])] )==0: continue
   J_year[ilat,ilon,ipft]=begy+np.nanargmin(J[ilat,ilon,ipft,:], axis=-1)
   if np.all(J[ilat,ilon,ipft,:]==J[ilat,ilon,ipft,0]):
    J_year[ilat,ilon,ipft]=endy

###Year to choose when nan
year_fill_J=Counter( J_year[J_year!=0].flatten()).most_common(1)[0][0]
J_year[J_year==0]=year_fill_J
np.save(homedir+"/"+"J_year",J_year)

#List the combinaison of the dimensions
list_dim=[]
name_restart={"SRF":"sechiba","SBG":"stomate"}
for rr in name_restart.keys():
 #Load a restart file as a template
 file_restart_0=dirout+"/../../../"+rr+"/Restart/"+exp2+"_"+str(begy)+"1231_"+name_restart[rr]+"_rest.nc"
 restart=xr.open_dataset(file_restart_0,decode_times=False,decode_cf=False)
 for var in restart.keys():
   list_dim.append(restart[var].dims)
list_dim=[xx for xx in list_dim if ('x' in xx)|('y' in xx)|('z_a' in xx)|('l_e' in xx)]
list_dim=np.unique(np.asarray(list_dim,dtype=object))

list_dim=list([list(xx) for xx in list_dim])
for ix,xx in enumerate(list_dim):
  for il,ll in enumerate(xx):
    if ll not in ["x","y","z_a","l_e"]:
     list_dim[ix][il]="ZZ"
list_dim=np.unique(np.array(list_dim,dtype=object))
print("Number of different dimensions: ",len(list_dim))

for ix,xx in enumerate(list_dim):
  for il,ll in enumerate(xx):
    if "l_e" in ll:
     list_dim[ix][il]="z_a"
list_dim=np.unique(list_dim)
print("Number of different dimensions 2: ",len(list_dim))

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
   #Open the restart file for the spinup year yy 
   f = Dataset(file_restart)
   ncf=f.variables[var][:]
   dimensions=list(f.variables[var].dimensions)
   f.close()
   for ix,xx in enumerate(dimensions):
    if xx not in ['x','y','l_e','z_a']:
     dimensions[ix]="ZZ"
    elif xx=="l_e":
     dimensions[ix]="z_a"
   dim_var=np.shape(ncf)
   if (dim_var==(1,1))|(len(dim_var)==1): continue
   if   dimensions == ('y','x'):
    restart[var].values[id_lat,id_lon]=np.copy(ncf[id_lat,id_lon])
   elif dimensions == ('ZZ','y','x'):
    restart[var].values[:,id_lat,id_lon]=np.copy(ncf[:,id_lat,id_lon])
   elif dimensions == ('ZZ','z_a','y','x'):
    restart[var].values[:,id_pft,id_lat,id_lon]=np.copy(ncf[:,id_pft,id_lat,id_lon])
   elif dimensions == ('ZZ','ZZ','y','x'): 
    restart[var].values[:,:,id_lat,id_lon]=np.copy(ncf[:,:,id_lat,id_lon])
   elif dimensions == ('ZZ','z_a','ZZ','y','x'):
    restart[var].values[:,id_pft,:,id_lat,id_lon]=np.copy(ncf[:,id_pft,:,id_lat,id_lon])
   elif dimensions == ('ZZ','ZZ','z_a','y','x'):
    restart[var].values[:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,id_pft,id_lat,id_lon])
   elif dimensions == ('ZZ','ZZ','ZZ','y','x'):
    restart[var].values[:,:,:,id_lat,id_lon]=np.copy(ncf[:,:,:,id_lat,id_lon])
   elif dimensions == ('ZZ','ZZ','ZZ','z_a','y','x'):
    restart[var].values[:,:,:,id_pft,id_lat,id_lon]=np.copy(ncf[:,:,:,id_pft,id_lat,id_lon])
   elif dimensions == ('ZZ','ZZ','z_a','ZZ','y','x'):
    restart[var].values[:,:,id_pft,:,id_lat,id_lon]=np.copy(ncf[:,:,id_pft,:,id_lat,id_lon])
   elif dimensions == ('ZZ','ZZ','ZZ','ZZ','y','x'):
    restart[var].values[:,:,:,:,id_lat,id_lon]=np.copy(ncf[:,:,:,:,id_lat,id_lon])

 os.system("rm -f "+homedir+"/"+name_restart[rr]+"_"+str(nb_ac)+"ac_f.nc")
 os.system("rm -f "+homedir+"/"+name_restart[rr]+"_"+str(nb_ac)+"ac_i.nc")
 restart.to_netcdf(homedir+"/"+name_restart[rr] +"_"+str(nb_ac)+"ac_f.nc")
 os.system("cp -r "+file_restart+ " "+homedir+"/"+name_restart[rr]+"_"+str(nb_ac)+"ac_i.nc")

#CREATION OF THE OUTPUT FILE------------------------------------------------------------------------------------------

print('Creation of the output files')
#Load an ORCHIDEE outputt file as a template
file_output_0=dirout+exp2+"_"+period+"_1M_stomate_history.nc"
restart=xr.open_dataset(file_output_0,decode_times=False,decode_cf=False)

#List of the dimensions
list_dim=[]
for var in restart.keys():
 list_dim.append(list(restart[var].dims))

for ix,xx in enumerate(list_dim):
  for il,ll in enumerate(xx):
    if ll not in ["veget","lon","lat"]:
     list_dim[ix][il]="ZZ"
list_dim=np.unique(np.array(list_dim, dtype=object))
list_dim=[xx for xx in list_dim if ("veget" in xx)|("lon" in xx)|("lat" in xx)]
#print(list_dim)

list_vars=["HEIGHT_DOM","DIA_DOM","VEGET_MAX","AGE","AGE_STAND"]
for yy in range(begy,endy+1):
  period=str(yy)+"0101_"+str(yy)+"1231"
  file_output=dirout+exp2+"_"+period+"_1M_stomate_history.nc"
  id_yr=np.where(J_year==yy)
  id_lon=id_yr[1] 
  id_lat=id_yr[0]
  id_pft=id_yr[2]
  if len(id_yr)==0: continue
  for var in list_vars:
   f = Dataset(file_output)
   ncf=f.variables[var][:]
   dimensions=list(f.variables[var].dimensions)
   f.close()
   for ix,xx in enumerate(dimensions):
    if xx not in ['lon','lat','veget']:
     dimensions[ix]="ZZ"
   dim_var=np.shape(ncf)
   if (dim_var==(1,1))|(len(dim_var)==1): continue
   if   dimensions == ('lat','lon'):
    restart[var].values[id_lat,id_lon]=np.copy(ncf[id_lat,id_lon])
   elif dimensions == ('ZZ','lat','lon'):
    restart[var].values[:,id_lat,id_lon]=np.copy(ncf[:,id_lat,id_lon])
   elif dimensions == ('ZZ','veget','lat','lon'):
    restart[var].values[:,id_pft,id_lat,id_lon]=np.copy(ncf[:,id_pft,id_lat,id_lon])
   elif dimensions == ('ZZ','ZZ','lat','lon'):
    restart[var].values[:,:,id_lat,id_lon]=np.copy(ncf[:,:,id_lat,id_lon])
os.system("rm -f "+homedir+"/output_"+str(nb_ac)+"ac_f.nc")
os.system("rm -f "+homedir+"/output_"+str(nb_ac)+"ac_i.nc")
    
restart.to_netcdf(homedir+"/output_"+str(nb_ac)+"ac_f.nc")
os.system("cp -r "+file_output+ " "+homedir+"/output_"+str(nb_ac)+"ac_i.nc")

