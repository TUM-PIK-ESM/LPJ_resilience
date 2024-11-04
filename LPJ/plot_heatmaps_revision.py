#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 19:11:04 2022

@author: bathiany
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import assign_properties
#import copy
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

### input
dpi=500


explist=('sba0082','sba0084', 'sba0087')
explist=('sba0084',)

varlist=('vegc',) #'agb','agb_tree')
#typelist=('ACyL10',)
#typelist=('SDy',)
typelist=('timmean',)


#explist=('sba0084',)
#varlist=('vegc','agb','agb_tree')
#typelist=('ACyL10',)



#explist=('sba0042','sba0044')
#varlist=('vegc','fpc','fpc_woody')
#typelist=('timmean',)

# 4 lags
#explist=('sba0082','sba0084',)
#varlist=('vegc',)
#typelist=('ACyL01','ACyL05','ACyL10','ACyL30',)
##typelist=('SDy',)


#varlist=('pft_cleafperm2','pft_crootperm2','pft_csapwperm2','pft_chawoperm2')
#varlist=('pft_frac_csapw',) #'pft_frac_chawo',
#typelist=('timmean',)
levellist_orig=[0,1]
#levellist_orig=[0,1,2,3,4,5,6,7,8,9]


fontsize_ticks=23
fontsize_title=12
fontsize_labels=30




levellist_orig=[0]
#levellist_orig=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]


for type in typelist:    
    if type == '':
        typestring=''
    else:
        typestring='_'+type

    for var in varlist:
        #print(var)
        varname, unit, minval, maxval, cmap_var=assign_properties.assign_varname(var,type)
    
        for exp in explist:

            #Tvec=assign_properties.assign_Tvec(exp) #[0]
            MAPvec=assign_properties.assign_Pvec(exp)            
            PSvec=assign_properties.assign_PSvec(exp) #[0]
            PSfactvec=assign_properties.assign_PSfactvec(exp) #[0]

            file = nc.Dataset(exp+'/'+var+typestring+'.nc')
            
            data = file[varname]

            ### take out zeros (not def)
            #data2=copy.copy(data) # fails
            #data2[data2==0]='NaN'
            
            dims=np.shape(np.shape(data))[0]
            
            if dims==4:
                levellist=levellist_orig
            else: 
                levellist=[0]
                        
            for level in levellist:

                plt.figure()
                #plt.rcParams.update({'font.size': 30})
                fig, ax = plt.subplots(figsize=(14,10))
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
                    else:    
                        Panel=plt.imshow(data[0,:,:], cmap=cmap_var, vmin=minval, vmax=maxval)

                
                ax.invert_yaxis()
                ax.set_yticks(np.arange(0,len(PSvec)))
                ax.set_xticks(np.arange(0,len(MAPvec)))
                #if dims==4:
                #    plt.title(exp+ ', ' + var+', index '+str(level)+', ' + type,  fontsize=fontsize_title)
                #else:
                #    plt.title(exp+ ', ' + var +', ' + type, fontsize=fontsize_title)
                    
                ax.set_xticklabels(MAPvec, fontsize=fontsize_ticks)
                ax.set_yticklabels(PSfactvec, fontsize=fontsize_ticks)
                cbar = plt.colorbar(Panel)
                if type.find('AC') != -1:
                    unit='Autocorrelation'
                if type == "timstd" or type== 'SDy':
                    unit='Standard deviation [gC/m²]'
                cbar.set_label(unit, rotation=90, fontsize=fontsize_labels)
                cbar.ax.tick_params(labelsize=fontsize_ticks)
                ax.set_ylabel("Precipitation amplitude factor [mm/yr]", fontsize=fontsize_labels)
                ax.set_xlabel("Precipitation mean [mm/yr]", fontsize=fontsize_labels)

                if dims==4:
                    figname=exp+'_'+var+typestring+'_tile'+str(level)+'_heatmap.png'
                else:                    
                    figname=exp+'_'+var+typestring+'_heatmap.png'
                
                if int(exp[5:7])>55:
                    ax.tick_params(axis='x', rotation=45)


                plt.savefig(figname, dpi=dpi)
                
                plt.show()
            
            del data
        
        
    