#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 18:36:59 2023

@author: bathiany
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import assign_properties
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

## plots a property versus precipitation from idealised LPJmL simulation


#explist=('sba0082','sba0084',)
explist=('sba0084',)

#typelist=('timmean','ACyL01','ACyL05','ACyL10',)

## Fig a, b (fpc and C)
typelist=('timmean',)
#varlist=('fpc_woody','vegc')

varlist=('pft_cleaf_frac','pft_croot_frac','pft_csapw_frac','pft_chawo_frac','pft_SH_frac')


## SD
#typelist=('SDy',)
#varlist=('vegc',)

### AC
#typelist=('ACyL10',)
#varlist=('vegc',)


## precip is on x axis
Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]
Pvec=np.array(Pvec)                 # in mm/year

Pvec=Pvec[1:20] # remove first value

dpi=500

fontsize_ticks=34   #23
fontsize_title=12
fontsize_labels=36  #30
fontsize_legend=34


for type in typelist:

    for var in varlist:
        varname, unit, minval, maxval, cmap_var=assign_properties.assign_varname(var,type)

        if type.find('AC') != -1:
            minval=0
            maxval=1
        elif type.find('decorr') != -1:
            minval=0
            maxval=110     
    
        for exp in explist:
         
            file = nc.Dataset(exp+'/'+var+'_'+type+'.nc')
            
            print(exp+'/'+var+'_'+type+'.nc')
    
            data = file[varname]
            print(np.shape(data))
            dims=np.shape(np.shape(data))[0]
            
            if dims==4:
                data_slice=data[0,0,0,1:20]
            else:
                data_slice=data[0,0,1:20]

            plt.figure()
            fig, ax = plt.subplots(figsize=(18,10))
            plt.plot(Pvec,data_slice, linewidth=5, color='k')
            
            varname_label=varname
            if var=='vegc':
                varname_label='C' + ' ['+unit+']'
            elif var=='fpc' or var=='fpc_woody':
                varname_label='FPC' + ' ['+unit+']'
            if type.find('AC') != -1:
                varname_label='autocorrelation'

            ax.set_ylabel(varname_label, fontsize=fontsize_labels)
            ax.set_xlabel('Precipitation [mm/yr]', fontsize=fontsize_labels)
    
            if var=='vegc' and type=='timmean':
                maxval=25000
            if var=='vegc' and type=='SDy':
                maxval=800

    
            plt.ylim([minval, maxval])
            plt.xlim([0, 2500])

            ax.tick_params(labelsize=fontsize_ticks)
    
            figname=exp+'_'+var+'_'+type+'_vs_Prec_crosssection.png'

            plt.savefig(figname, dpi=dpi)
    
            plt.show()








