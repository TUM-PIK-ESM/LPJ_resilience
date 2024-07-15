#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 09:19:05 2024

@author: bathiany
"""

import netCDF4 as nc
#import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from matplotlib.colors import LogNorm
import assign_properties


dpi=1000
fontsize_ticks=18
fontsize_title=22
fontsize_labels=18


# local:
#sys.path.insert(1, '/home/bathiany/Projects/Vegetation_resilience_indicators')

explist=('sba0080',)
varlist=('fpc',)
typelist=('timmean',)

levellist=[0,1]


for type in typelist:    
    if type == '':
        typestring=''
    else:
        typestring='_'+type

    for var in varlist:
        #print(var)
        varname, unit, minval, maxval, cmap_var=assign_properties.assign_varname(var,type)
    
        for exp in explist:
            
            file = nc.Dataset(exp+'/'+var+typestring+'.nc')
            
            data = file[varname]
            
##          print(np.shape(data))
            lons = file.variables['longitude'][:]
            lats = file.variables['latitude'][:]
                        
            dims=np.shape(np.shape(data))[0]
            #print(dims)
            if dims==4:
                print("4 dimensions detected")

            data2=np.squeeze(data)
            ## mask 0
            data2[data2==0]=np.nan
            # misval:
            data2[data2<-10**31]=np.nan
            
            for level in levellist:    
                
                #### plot with basemap
                plt.figure(dpi=dpi)
                plt.rcParams.update({'font.size': 30})
                fig, ax = plt.subplots(figsize=(6,7))
    
                print(np.shape(data2))
                #print(np.shape(data2[0,level,:,:]))
                
                # Find the range of longitudes to center the plot around the zero meridian
                lon_range = max(lons) - min(lons)
                center_lon = (max(lons) + min(lons)) / 2.0
    
                # Adjust the longitudes to center the plot around the zero meridian
                lons_centered = lons - center_lon
                lons_centered[lons_centered < -180] += 360
    
                m = Basemap(projection='cyl', llcrnrlon=center_lon - lon_range / 2, llcrnrlat=min(lats), urcrnrlon=center_lon + lon_range / 2, urcrnrlat=max(lats), resolution='l')
                #m = Basemap(projection='cyl', llcrnrlon=min(lons_centered), llcrnrlat=min(lats), urcrnrlon=max(lons_centered), urcrnrlat=max(lats), resolution='l') 
                x, y = m(lons, lats)
                #x, y = m(lons_centered, lats)                    
                           
                m.drawcoastlines()
                m.drawcountries()
                if type.find('AC') != -1:
                    print('AC detected')
                    vmin=0
                    vmax=1
                    data_clipped = np.clip(data2[level,:,:], vmin, vmax)
                    m.pcolormesh(x, y, data_clipped, cmap='viridis', shading='auto', vmin=vmin, vmax=vmax)   
                    cbar = plt.colorbar(label='', orientation='horizontal', pad=0.05, shrink=0.75)
                    cbar.ax.tick_params(labelsize=fontsize_ticks)  # Change font size here
                    cbar.set_label('Autocorrelation', fontsize=fontsize_labels)
                
                else:
                    vmin=0.2
                    vmax=0.8
                    m.pcolormesh(x, y, data2[level,:,:], cmap='viridis', shading='auto', vmin=vmin, vmax=vmax)
                    cbar = plt.colorbar(label='', orientation='horizontal', pad=0.05, shrink=0.75)
                    cbar.ax.tick_params(labelsize=fontsize_ticks)  # Change font size here

                    cbar.set_label(varname, fontsize=fontsize_labels)
                            
                #plt.show()
    
                figname=exp+'_'+var+typestring+'_tile'+str(level)+'_geographical_map.png'
    
                plt.savefig(figname, dpi=dpi)
    
                
                plt.show()
    #   
    #                plt.tight_layout()
