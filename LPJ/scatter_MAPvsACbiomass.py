#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 18:40:11 2022

@author: bathiany
"""

import matplotlib.pyplot as plt
import netCDF4 as nc
from os.path import exists
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

dpi=1000
size=2

fontsize_ticks=18
fontsize_title=22
fontsize_labels=24

############

exp='sba0080'

case_dict = {
 '1': ('Prec','Prec','timmean','mm/yr',3600*24*365, 0, 4000, 'agb_tree','agb_tree', 'ACyL10', '', 1, 0, 1), 
 '2': ('Tair','Tair','timmean','K',1, 285, 302, 'agb_tree','agb_tree', 'ACyL10', '', 1, 0, 1),
 '3': ('LWdown','LWdown','timmean','W/m2',1, 200, 440, 'agb_tree','agb_tree', 'ACyL10', '', 1, 0, 1),
 '4': ('SWdown','SWdown','timmean','W/m2',1, 160, 270, 'agb_tree','agb_tree', 'ACyL10', '', 1, 0, 1),
 '5': ('Prec','Prec','anncycamp','mm/yr',3600*24*365, 0, 3000, 'agb_tree','agb_tree', 'ACyL10', '', 1, 0, 1), 
 '6': ('Prec','Prec','daystd','mm/yr',3600*24*365, 0, 3000, 'agb_tree','agb_tree', 'ACyL10', '', 1, 0, 1), 
 '7': ('Prec','Prec','yearstd','mm/yr',3600*24*365, 0, 250, 'agb_tree','agb_tree', 'ACyL10', '', 1, 0, 1), 

 }

case=1
axesranges_on=1


var1=case_dict[str(case)][0]
varname1=case_dict[str(case)][1]
type1=case_dict[str(case)][2]
unit1=case_dict[str(case)][3]
scale1=case_dict[str(case)][4]
xmin=case_dict[str(case)][5]
xmax=case_dict[str(case)][6]
var2=case_dict[str(case)][7]
varname2=case_dict[str(case)][8]
type2=case_dict[str(case)][9]
unit2=case_dict[str(case)][10]
scale2=case_dict[str(case)][11]
ymin=case_dict[str(case)][12]
ymax=case_dict[str(case)][13]


#########

# in
file1 = '/home/bathiany/tipes/input3/'+var1+'_WFDE5_CRU_v1.0_daily_resampled_1000yrs_mesosouthamerica_fixed_'+type1+'.nc'
file2 = exp+'/'+var2+'_'+type2+'.nc'

# out
file_fig=exp+'_'+var1+'_'+type1+'_'+var2+'_'+type2+'_scatter.png'
  

if exists(file1) and exists(file2):
    
    data1 = nc.Dataset(file1)
    data2 = nc.Dataset(file2)
    vardata1 = data1[varname1][:,:]    
    vardata2 = data2[varname2][:,:] 

    vardata1[vardata1==1]='NaN'
    vardata1[vardata1==-1]='NaN'
    vardata2[vardata2==1]='NaN'
    vardata2[vardata2==-1]='NaN'
    
    plt.figure(dpi=dpi)
    plt.rcParams.update({'font.size': 30})
    fig, ax = plt.subplots(figsize=(10,7))
            
    plt.scatter(vardata1*scale1,vardata2*scale2, s=size, c='k')
    if axesranges_on==1:
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

    if type2.find('AC') != -1:
        plt.ylabel('Autocorrelation', fontsize=fontsize_labels)
    else:
        plt.ylabel(var2 + ' [' + unit2 + '], '+type2, fontsize=fontsize_labels)

    if var1=='Prec':
        plt.xlabel('Precipitation [mm/yr]', fontsize=fontsize_labels)
    elif var1=='Tair':
        plt.xlabel('Temperature [K]', fontsize=fontsize_labels)

    plt.savefig(file_fig)
    plt.show()
    

