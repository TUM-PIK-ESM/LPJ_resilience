# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats, scipy.signal, scipy.optimize

### estimates recovery rate, based on Taylor Smith's code; 
## see Smith, Traxl, Boers; Nature Clim Change, 2022
## modified by Sebastian Bathiany, 24.11.24.

p0 = [-5000,-1/50]  ## initial values for parameters a and b in f(t)=a*exp(b*t)
bounds = ([-100000, -1], [100000, 0])  # bounds of params a,b  ;  a < 0, b < 0 for recovery from below


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

