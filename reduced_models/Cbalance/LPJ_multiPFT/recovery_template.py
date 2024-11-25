#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 19:46:43 2024

@author: bathiany
"""

import netCDF4 as nc
import sys
import numpy as np
import subprocess

### edit
exp='experiment_id'
var='var_id'
yr_perturb=yr_perturb_id

margin=1   # stop taking data into account if closer to 0 than this margin
## in terms of state? or in terms of time (datapoints) ?

def write_data(data, nc_filename, varname):

    nc_file = nc.Dataset(nc_filename, mode="a")

    var_nc = nc_file.variables[varname]
    
    # Write data into the variable
    var_nc[0,:, :] = data
    
    # Close the NetCDF file
    nc_file.close()
    

Vegpath='/home/bathiany/Projects/Vegetation_resilience_indicators/'
LPJmultipath=Vegpath+'/reduced_models/Cbalance/LPJ_multiPFT/'
LPJpath=Vegpath + 'LPJ/'
functionspath=Vegpath + 'reduced_models/Cbalance/functions/'

sys.path.insert(1, functionspath)
import LPJ_functions

sys.path.insert(1, Vegpath)
import empirical_recovery_SB

sys.path.insert(1, LPJpath)
import assign_properties

## create templatefile
templatefile=LPJmultipath + exp + '/'+var+'.nc'

#command = ['cdo mulc,0 -sellevel,1/'+str(NPFT)+' -selyear,'+str(y_ini)+'/'+str(y_fin) + ' ' + templatefile + ' ' + outfile]

outfile_recovtime=LPJmultipath + exp + '/'+var+'_recoverytime.nc'
outfile_recovrate=LPJmultipath + exp + '/'+var+'_recoveryrate.nc'
outfile_rsq=LPJmultipath + exp + '/'+var+'_rsq.nc'
command1 = ['cdo mulc,0 -seltimestep,1 ' + templatefile + ' ' + outfile_recovtime]
command2 = ['cdo mulc,0 -seltimestep,1 ' + templatefile + ' ' + outfile_recovrate]
command3 = ['cdo mulc,0 -seltimestep,1 ' + templatefile + ' ' + outfile_rsq]

result1 = subprocess.run(command1, shell=True, capture_output=True, text=True)
result2 = subprocess.run(command2, shell=True, capture_output=True, text=True)
result3 = subprocess.run(command3, shell=True, capture_output=True, text=True)

####

varname=assign_properties.assign_varname(var,'')[0]

file = nc.Dataset(LPJmultipath+exp+'/vegc_anom2stat.nc')
var_all_recovery = file[varname]


dims=var_all_recovery.shape
ndims=np.size(dims)
if ndims==3:
    var_recovery=LPJ_functions.unmask(var_all_recovery[:,:,:])

elif ndims==4:
    var_recovery=LPJ_functions.unmask(var_all_recovery[:,:,:,:])

rate_matrix=np.zeros((dims[1],dims[2]))
time_matrix=np.zeros((dims[1],dims[2]))
rsq_matrix=np.zeros((dims[1],dims[2]))

## open outfile for writing
for lat in range(0,dims[1]):
    for lon in range(0,dims[2]):
        ## for level in ...
        a, b, rsq = empirical_recovery_SB.fit_to_transition_SB(var_recovery[yr_perturb+1:,lat,lon])
        time_matrix[lat,lon]=-1/b
        rate_matrix[lat,lon]=-b
        rsq_matrix[lat,lon]=rsq

write_data(rate_matrix, outfile_recovrate, varname)
write_data(time_matrix, outfile_recovtime, varname)
write_data(rsq_matrix, outfile_rsq, varname)
