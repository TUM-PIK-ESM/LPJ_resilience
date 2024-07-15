#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:25:30 2023

@author: bathiany
"""
### function for computing AC
##https://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
def ACF_fun(x,lagmax):
    import numpy as np
    import math
    mean=x.mean()
    var=np.var(x)
    xp=x-mean
    if math.isnan(mean) or math.isnan(var):
        #print('fail')
        corr=np.empty((lagmax))
        corr[:] = np.nan
    else:
        corr=np.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)
    
    return corr[0:lagmax]


## running window
def ts_separate(ts, bw):
    # applies smoothing with certain bandwidth bw
    # then removes the smoothed low-freq component

    import numpy as np
    import pandas as pd

    df = pd.DataFrame(data=ts)

    # https://medium.com/@jodancker/a-brief-introduction-to-time-series-smoothing-4f7ed61f78e1
    ts_smoothed = df.rolling(window=bw, win_type="gaussian", center=True).mean(std=bw/2.3548)
    
    ts_smoothed=ts_smoothed.values[:,0]   # convert back to array

    ts_masked=ts.copy()
    ts_masked[np.isnan(ts_smoothed)] = ts_smoothed[np.isnan(ts_smoothed)]

    ts_anomaly=ts_masked-ts_smoothed

    ts_anomaly = ts_anomaly[~np.isnan(ts_anomaly)]

    return ts_anomaly, ts_smoothed