#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 18:40:11 2022

@author: bathiany
"""

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

dpi=1000
size=2

fontsize_ticks=18
fontsize_title=22
fontsize_labels=24

############

exp='sba0080'



case=2
axesranges_on=1


#########
# in
Pfile = '/home/bathiany/tipes/input3/Prec_WFDE5_CRU_v1.0_daily_resampled_1000yrs_mesosouthamerica_fixed_timmean.nc'
Tfile = '/home/bathiany/tipes/input3/Tair_WFDE5_CRU_v1.0_daily_resampled_1000yrs_mesosouthamerica_fixed_timmean.nc'
Cfile = exp+'/agb_tree_ACyL10.nc'


Pdata = nc.Dataset(Pfile)
Tdata = nc.Dataset(Tfile)
Cdata = nc.Dataset(Cfile)

Pvardata = Pdata['Prec'][:,:]
Tvardata = Tdata['Tair'][:,:]
Cvardata = Cdata['agb_tree'][:,:]

Pvardata[Pvardata==1]='NaN'
Pvardata[Pvardata==-1]='NaN'
Tvardata[Tvardata==1]='NaN'
Tvardata[Tvardata==-1]='NaN'
Cvardata[Cvardata==1]='NaN'
Cvardata[Cvardata==-1]='NaN'


Pvardata=Pvardata*3600*24*365
Tvardata=Tvardata-273.15



plt.figure(dpi=dpi)
plt.rcParams.update({'font.size': 30})
fig, ax = plt.subplots(figsize=(10,7))

if case==1:
    file_fig=exp+'_Prec_AGBtree_scatter_colored3d.png'
    sc=plt.scatter(Pvardata,Cvardata, s=size, c=Tvardata, cmap='viridis', vmin=0, vmax=30)
    plt.xlim(0, 4000)

    plt.xlabel('Precipitation [mm/yr]', fontsize=fontsize_labels)

    # Add colorbar with label
    cbar = plt.colorbar(sc)
    cbar.set_label('Temperature [°C]', fontsize=fontsize_labels)  # Adjust label as needed
    
    # Customize the colorbar ticks (optional)
    cbar.set_ticks(np.linspace(0, 30, 5))  # Example: 4 ticks


elif case==2:
    file_fig=exp+'_Tair_AGBtree_scatter_colored3d.png'
    sc=plt.scatter(Tvardata,Cvardata, s=size, c=Pvardata, cmap='viridis', vmin=0, vmax=4000)
    plt.xlim(12, 30)
    
    plt.xlabel('Temperature [°C]', fontsize=fontsize_labels)

    # Add colorbar with label
    cbar = plt.colorbar(sc)
    cbar.set_label('Precipitation [mm/yr]', fontsize=fontsize_labels)  # Adjust label as needed
    
    # Customize the colorbar ticks (optional)
    cbar.set_ticks(np.linspace(0, 4000, 5))  # Example: 4 ticks

plt.ylim(0, 1)

plt.ylabel('Autocorrelation', fontsize=fontsize_labels)

plt.savefig(file_fig)
plt.show()


