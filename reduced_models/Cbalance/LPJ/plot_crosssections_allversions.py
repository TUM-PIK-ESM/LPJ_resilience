#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 18:36:59 2023

@author: bathiany
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(1, '/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ')
#import netCDF4 as nc
import assign_properties
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

## requires reading npy files generated by Cbalance_model_sensitivity_revised.py

exp="sba0084"      # exp from which data is read in
var='vegc'
#versions=("A1N1M1E1A1S1","A1N1M0E1A1S1","A1N3M0E1A1S1","A2N1M0E1A1S1")
#labels=("LPJ-CNi","LPJ-CN","LPJ-C","LPJ-N")

versions=("A1N1M0E1A1S1","A1N3M0E1A1S1","A2N1M0E1A1S1")
labels=("LPJ-CN","LPJ-C","LPJ-N")
type='AC10_Cperm2'
#type='AC10_N'




## The versions for overview plot:
# orig reduced model = A1N1M1
# LPJ-CN = A1N1M0
# LPJ-C = A1N3M0
# LPJ-N = A2N1M0


parameter="varfrac"
paramind=5 #default


## precip is on x axis
Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]
Pvec=np.array(Pvec)                 # in mm/year

#Pvec=Pvec[1:20] # remove first value
cellind_ini=1   # remove P=0, starting from 1 (200 mm/yr, i.e. cellind=1)
cellind_fin=20
Pvec=Pvec[cellind_ini:cellind_fin] 

Nprec=cellind_fin-cellind_ini


## the sensitivity script only simulates a reduced range of Pvec, not the full range.
#Nprec=10   # number of grid cells, starting from 1 (200 mm/yr, i.e. cellind=1)
#Pvec=Pvec[1:Nprec+1] # remove P=0, and cut off values above 2000mm



dpi=500

fontsize_ticks=34 #23
fontsize_title=12
fontsize_labels=36 #30
fontsize_legend=34



varname, unit, minval, maxval, cmap_var=assign_properties.assign_varname(var,type)

if type.find('AC') != -1:
    minval=-0.1
    maxval=1
elif type.find('decorr') != -1:
    minval=0
    maxval=110

plt.figure()
fig, ax = plt.subplots(figsize=(18,10))



vind=0
for version in versions:
    filename=exp + '_' + version + '_'+parameter+'_vec_'+type+'_diag.npy'
    if os.path.isfile(filename):
        data=np.load(filename)
        data=data[paramind,cellind_ini:cellind_fin]
        #print(data)
        #print(filename)
    #print(np.shape(data))
    
    plt.plot(Pvec,data, label=labels[vind], linewidth=5)
    vind=vind+1

varname_label=varname
if var=='vegc':
    varname_label='C' + ' ['+unit+']'
elif var=='fpc' or var=='fpc_woody':
    varname_label='FPC' + ' ['+unit+']'
if type.find('AC') != -1:
    varname_label='autocorrelation'

ax.set_ylabel(varname_label, fontsize=fontsize_labels)
ax.set_xlabel('Precipitation [mm/yr]', fontsize=fontsize_labels)

plt.xlim([0, 2500])

plt.ylim([minval, maxval])
if type.find('AC') != -1:
    ax.legend(fontsize=fontsize_legend, loc='upper left') # upper right
else:
    ax.legend(fontsize=fontsize_legend, loc='lower right')

#       ax.set_xticklabels(Pvec, fontsize=fontsize_ticks)
#       ax.set_yticklabels(fontsize=fontsize_ticks)
ax.tick_params(labelsize=fontsize_ticks)


figname=var+'_'+type+'_vs_Prec_crosssection_3versions.png'
#figname=var+'_'+type+'_vs_Prec_crosssection_4versions.png'


plt.savefig(figname, dpi=dpi)

plt.show()


