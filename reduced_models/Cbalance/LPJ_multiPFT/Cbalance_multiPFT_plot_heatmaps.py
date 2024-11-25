#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: bathiany
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import assign_properties
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')


dpi=500

levellist_orig=[0,1]

fontsize_ticks=18
fontsize_title=20
fontsize_labels=18

varlist=('vegc',)

explist=('sba0056_P0000',)
typelist=('ACyL10',)

#explist=('sba0056_P0010', 'sba0056_P0013', 'sba0056_P0014', 'sba0056_P0015', 
#         'sba0056_P0020', 'sba0056_P0023', 'sba0056_P0024', 'sba0056_P0025', 
#         'sba0056_P0030', 'sba0056_P0033', 'sba0056_P0034', 'sba0056_P0035')

explist=('sba0056_P0024','sba0056_P0034')
typelist=('recoverytime','recoveryrate')
#typelist=('recoverytime',)


levellist_orig=[0]


for var in varlist:
    #print(var)
    varname, unit, minval, maxval, cmap_var=assign_properties.assign_varname(var,type)

    for exp in explist:

        for type in typelist:    
            if type == '':
                typestring=''
            else:
                typestring='_'+type


            Tvec=assign_properties.assign_Tvec(exp) #[0]
            Pvec=assign_properties.assign_Pvec(exp)            
                      
            file = nc.Dataset(exp+'/'+var+typestring+'.nc')
            
            data = file[varname]
            
            dims=np.shape(np.shape(data))[0]
            

            if exp=='sba0035' or exp=='sba0036' or exp=='sba0044':
                levellist=[0]
            elif dims==4:
                levellist=levellist_orig
            else: 
                levellist=[0]
            
            
            for level in levellist:
                
                plt.figure()
                #plt.rcParams.update({'font.size': 30})
                #fig, ax = plt.subplots(figsize=(14,10))
                fig, ax = plt.subplots(figsize=(10,7))

                #fig, ax = plt.subplots()
                
                
                if dims==4:                    
                    if type.find('AC') != -1:
                        Panel=plt.imshow(data[0,level,:,:], cmap=cmap_var, vmin=0, vmax=1)
                    elif type.find('decorr') != -1:
                        Panel=plt.imshow(data[0,level,:,:], cmap=cmap_var, vmin=0, vmax=150)
                    else:                      
                        Panel=plt.imshow(data[0,level,:,:], cmap=cmap_var, vmin=minval, vmax=maxval)

                else:
                    if type.find('AC') != -1:
                        Panel=plt.imshow(data[0,:,:], cmap=cmap_var, vmin=0, vmax=1)
                    elif type.find('decorr') != -1:
                        Panel=plt.imshow(data[0,:,:], cmap=cmap_var, vmin=0, vmax=150)
                    elif type == 'corr':
                        Panel=plt.imshow(data[0,:,:], cmap=cmap_var, vmin=-1, vmax=1)
                    elif type == "recoverytime":
                        Panel=plt.imshow(data[0,:,:], cmap=cmap_var, vmin=0, vmax=130)
                    elif type == "recoveryrate":
                        Panel=plt.imshow(data[0,:,:], cmap=cmap_var, vmin=0, vmax=0.1)
                    else:    
                        Panel=plt.imshow(data[0,:,:], cmap=cmap_var, vmin=minval, vmax=maxval)

                
                ax.invert_yaxis()
                ax.set_xticks(np.arange(0,len(Tvec)))
                ax.set_yticks(np.arange(0,len(Pvec)))
                if dims==4:
                    plt.title(exp+ ', ' + var+', index '+str(level)+', ' + type,  fontsize=fontsize_title)
                else:
                    plt.title(exp+ ', ' + var +', ' + type, fontsize=fontsize_title)
                    
                ax.set_xticklabels(Tvec, fontsize=fontsize_ticks)
                ax.set_yticklabels(Pvec, fontsize=fontsize_ticks)
                cbar = plt.colorbar(Panel)
                if type.find('AC') != -1:
                    unit='Autocorrelation'
                elif type == "timstd":
                    unit='Standard deviation [gC/m²]'
                elif type == "recoverytime":
                    unit='years'      
                elif type == "recoveryrate":
                    unit='1/years' 
                cbar.set_label(unit, rotation=90, fontsize=fontsize_labels)
                cbar.ax.tick_params(labelsize=fontsize_ticks)
                ax.set_xlabel("Temperature [°C]", fontsize=fontsize_labels)
                ax.set_ylabel("Precipitation [mm/yr]", fontsize=fontsize_labels)

                if dims==4:
                    figname=exp+'_'+var+typestring+'_tile'+str(level)+'_heatmap.png'
                else:                    
                    figname=exp+'_'+var+typestring+'_heatmap.png'
                
                if int(exp[5:7])>55:
                    ax.tick_params(axis='x', rotation=-45)

                plt.savefig(figname, dpi=dpi)
                
                plt.show()
            
            del data
        
        
    