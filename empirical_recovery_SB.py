# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats, scipy.signal, scipy.optimize
import netCDF4 as nc
import subprocess
import sys

### estimates recovery rate, based on Taylor Smith's code; 
## see Smith, Traxl, Boers; Nature Clim Change, 2022
## modified by Sebastian Bathiany, 24.11.24.

p0 = [-5000,-1/50]  ## initial values for parameters a and b in f(t)=a*exp(b*t)
bounds = ([-100000, -1], [100000, 0])  # bounds of params a,b  ;  a < 0, b < 0 for recovery from below

        
def compute_recovery(var, datapath, yr_perturb):
    Vegpath='/home/bathiany/Projects/Vegetation_resilience_indicators/'
    LPJpath=Vegpath + 'LPJ/'
    functionspath=Vegpath + 'reduced_models/Cbalance/functions/'
    
    sys.path.insert(1, functionspath)
    import LPJ_functions
    
    sys.path.insert(1, Vegpath)
    
    sys.path.insert(1, LPJpath)
    import assign_properties
    
    yr_perturb = int(yr_perturb)  # Convert to int if it's a string
    
    ## create templatefile
    templatefile=datapath + '/'+var+'.nc'
    
    #command = ['cdo mulc,0 -sellevel,1/'+str(NPFT)+' -selyear,'+str(y_ini)+'/'+str(y_fin) + ' ' + templatefile + ' ' + outfile]
    
    outfile_recovtime=datapath + '/'+var+'_recoverytime.nc'
    outfile_recovrate=datapath + '/'+var+'_recoveryrate.nc'
    outfile_rsq=datapath + '/'+var+'_rsq.nc'
    command1 = ['cdo mulc,0 -seltimestep,1 ' + templatefile + ' ' + outfile_recovtime]
    command2 = ['cdo mulc,0 -seltimestep,1 ' + templatefile + ' ' + outfile_recovrate]
    command3 = ['cdo mulc,0 -seltimestep,1 ' + templatefile + ' ' + outfile_rsq]
    
    subprocess.run(command1, shell=True, capture_output=True, text=True)
    subprocess.run(command2, shell=True, capture_output=True, text=True)
    subprocess.run(command3, shell=True, capture_output=True, text=True)
    
    ####
    
    varname=assign_properties.assign_varname(var,'')[0]
    
    file = nc.Dataset(datapath+'/vegc_anom2stat.nc')
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
            a, b, rsq = fit_to_transition_SB(var_recovery[yr_perturb+1:,lat,lon])
            time_matrix[lat,lon]=-1/b
            rate_matrix[lat,lon]=-b
            rsq_matrix[lat,lon]=rsq
    
    write_data(rate_matrix, outfile_recovrate, varname)
    write_data(time_matrix, outfile_recovtime, varname)
    write_data(rsq_matrix, outfile_rsq, varname)




def write_data(data, nc_filename, varname):

    nc_file = nc.Dataset(nc_filename, mode="a")

    var_nc = nc_file.variables[varname]
    
    # Write data into the variable
    var_nc[0,:, :] = data
    
    # Close the NetCDF file
    nc_file.close()
    
    
## To determine up to which point we should use the data (time of full recovery)
## but need to work with a marging, to stop beforehand
def find_zero_crossing(data, margin_rel=0.05):
    # Ensure the input is a NumPy array
    data = np.array(data)
    
    # Here we assume the time series is already an anomaly to the steady state and approaches 0!
    # Hence, the absolute range of the data is equal to data[0]
    data_shifted=data-data[0]*margin_rel # shift toward zero by margin => detect 0-transition earlier
    
    # Check for zero crossing by comparing the sign of consecutive elements
    sign_changes = np.sign(data_shifted[:-1]) != np.sign(data_shifted[1:])
    
    # Find the index of the first sign change
    crossing_index = np.argmax(sign_changes)  # Returns the first True index
    
    # Check if a crossing was detected
    if sign_changes[crossing_index]:
        return crossing_index + 1  # Return the index where it crosses zero
    else:
        return None  # Return None if no crossing is found

def exp_fit(x, a, b):
    return a*np.exp(b*x)

# derivative of f w.r.t. a and w.r.t. b
def exp_jac(x, a, b):
    jac = np.array([np.exp(b*x), a*x*np.exp(b*x)]).T
    return jac

# computes goodness of fit of data to exp function
# uses exp_fit
def get_rsq(ser, popt):
    residuals = ser - exp_fit(np.arange(ser.shape[0]), *popt)
    ss_res = np.nansum(residuals**2)
    ss_tot = np.nansum((ser - np.nanmean(ser))**2)
    if ss_tot == 0:
        r_squared=np.nan
    else:
        r_squared = 1 - (ss_res / ss_tot)
    return r_squared

# Get the fit statistics
def fit_to_transition_SB(data):
    ## use only the part of the data that is still recovering:
    index_fin=find_zero_crossing(data)
    data=data[0:index_fin]
    
    try:
        trange = np.arange(data.shape[0])
        popt, _ = scipy.optimize.curve_fit(exp_fit, trange, data, p0=p0, jac=exp_jac, bounds=bounds)
        rsq = get_rsq(data, popt)
        fit_a = popt[0]
        fit_b = popt[1]
    except Exception as e:
        print(f"Error: {e}")
        rsq = np.nan
        fit_a = np.nan
        fit_b = np.nan

    # if r_squared was not possible to compute because ss_tot was 0:
    if rsq==np.nan:
        fit_a = np.nan
        fit_b = np.nan        

    # if estimated recovery time is longer than the state actually needed to recover,
    # something is very wrong => discard
    if np.abs(1/fit_b) > data.shape[0]:
        rsq = np.nan
        fit_a = np.nan
        fit_b = np.nan        

    return fit_a, fit_b, rsq


if __name__ == "__main__":
    if len(sys.argv) > 1:
        compute_recovery(sys.argv[1], sys.argv[2], sys.argv[3])