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

explist=('sba0082', 'sba0084',)
#explist=('sba0042',)  # all

typelist=('ACyL01','ACyL05','ACyL10','ACyL30')
#lags=('1 year', '5 years', '10 years', '30 years')
lags=('1 yr.', '5 yr.', '10 yr.', '30 yr.')

colors=('#66CCEE','#CCBB44','#EE6677','#AA3377')


var='vegc'



## precip is on x axis
Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]
Pvec=np.array(Pvec)                 # in mm/year

Pvec=Pvec[1:20] # remove first value

index=0 # index to select slice

dpi=500

fontsize_ticks=34 #23
fontsize_title=12
fontsize_labels=36 #30
fontsize_legend=30 #34


varname, unit, minval, maxval, cmap_var=assign_properties.assign_varname(var,type)

minval=-0.1
maxval=1   

for exp in explist:

    plt.figure()
    fig, ax = plt.subplots(figsize=(18,10))
     
    lagind=0
    for type in typelist:

        file = nc.Dataset(exp+'/'+var+'_'+type+'.nc')
        
        data = file[varname]
        
        data_slice=data[0,index,1:20]
                
        plt.plot(Pvec,data_slice, label='lag: ' + lags[lagind], linewidth=5, color=colors[lagind])
        lagind+=1
    
    varname_label='autocorrelation'

    ax.set_ylabel(varname_label, fontsize=fontsize_labels)
    ax.set_xlabel('Precipitation [mm/yr]', fontsize=fontsize_labels)

    plt.ylim([minval, maxval])
    plt.xlim([0, 2500])

    ax.legend(fontsize=fontsize_legend, loc='lower left')


#            ax.set_xticklabels(Pvec, fontsize=fontsize_ticks)
#            ax.set_yticklabels(fontsize=fontsize_ticks)
    ax.tick_params(labelsize=fontsize_ticks)

    figname=exp+'_'+var+'_AC_vs_Prec_crosssection_multilags.png'


    plt.savefig(figname, dpi=dpi)

    plt.show()








