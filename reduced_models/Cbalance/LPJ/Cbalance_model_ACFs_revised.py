#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:51:54 2023

@author: bathiany
"""

### forced by NPP (bm_inc) and soil moisture from LPJ
## includes allocation_tree module

## Warning:
# with allocation enabled, model does not give identical results to LPJ anymore,
# at least at dry grid cells.
# However, model stays stable, and correlation between LPJ and this implementation is around ~0.99
# The mean values can be a bit off ( 8-9% for chawo, rest smaller)

import numpy as np
#import matplotlib as mpl
#mpl.rcParams.update(mpl.rcParamsDefault)
import matplotlib.pyplot as plt

import netCDF4 as nc
#from scipy import stats
#import math
import random
#from sklearn.metrics import mean_squared_error
import sys
sys.path.insert(1, '/home/bathiany/Projects/Vegetation_resilience_indicators/reduced_models/Cbalance/functions')
from ACF import ACF_fun
import LPJ_functions
import par

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

LPJdatapath='/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/'


########## parameters
dpi=500
exp="sba0084"

lags=100

cellindlist=[10,2]
#cellindlist=[6,7,8,9,10,11]

# wet cell, 2000 mm, has cellind 10,
# dry cell, 400 mm, has cellind 2

index_PS=0

Nparams=7
PFT=1

# default
parameter="varfrac"
paramind=5

#parameter="varfrac"
#paramind=0

plot_AC_single=0
plot_AC_cellpairs=0         # wet vs dry
plot_AC_versions=0          # version comparison
plot_AC_all=1               # all


#### figures for the paper

#nametag="TEST"
#versions=("A1N1M1E1A1S1",)

## Fig ACF a
#nametag="shuffling"
#versions=("A1N1M1E1A1S0","A1N1M1E1A1S1")

#nametag="NoNdyn"
#versions=("A1N1M1E1A1S1", "A1N3M0E1A1S1")

#nametag="NoNdyn2"
#versions=("A1N1M1E1A1S1", "A1N3M1E1A1S1")
#
#nametag="NoNdyn3"
#versions=("A1N1M1E1A1S1", "A1N0M0E1A1S1")

#nametag="NoNdyn4"
#versions=("A1N1M1E1A1S1", "A1N3M0E1A1S1", "A1N3M1E1A1S1",)

#nametag="shuffling_and_NoNdyn"
#versions=("A1N1M1E1A1S0","A1N1M1E1A1S1","A1N3M0E1A1S1")


#### demo, PIK
#nametag="PIK1"
#versions=("A1N1M1E1A1S1",)
#
#nametag="PIK2"
#versions=("A1N1M1E1A1S1","A1N3M0E1A1S1")
#
##nametag="PIK3"
##versions=("A1N1M1E1A1S0",)
##
##nametag="PIK4"
##versions=("A1N1M1E1A1S0","A1N3M0E1A1S0")


#nametag="AC_from_N_to_C"
#versions=("A1N2M0E1A1S0","A1N2M0E1A1S1")  # S0 has the AC pattern, S1 not
#
#nametag="AC_from_N_to_C2"
#versions=("A1N2M1E1A1S0","A1N2M1E1A1S1")  # S0 has the AC pattern, S1 not
#
#nametag="AC_from_N_to_C3"
#versions=("A1N3M0E1A1S1","A1N2M0E1A1S1")   ## all equal (due to shuffling in A1N2)
#
#nametag="AC_from_N_to_C4"
#versions=("A1N3M0E1A1S0","A1N2M0E1A1S0")   # A1N2 has the AC pattern, N3 not
#
## => The pattern is not present for any N3 version (fixed N)
## It appears for N2, but only if I do not shuffle.
#
### total versions:
## Fig ACF b
#nametag="AC_from_N_to_C5"
#versions=( "A1N3M0E1A1S0", "A1N3M0E1A1S1","A1N2M0E1A1S0","A1N2M0E1A1S1")
# called CmNf_R0 CmNf_R1 CiNp_R0 CiNp_R1 in submitted version

## revised paper, Fig ACF1
#nametag="AC1"
#versions=("A1N1M1E1A1S1","A1N3M0E1A1S1",)
# LPJ-CN with interactive mort and LPJ-C

#nametag="AC1_CN_vs_C"
#versions=("A1N1M0E1A1S1","A1N3M0E1A1S1",)

nametag="AC1_shuffle_and_mort"
versions=("A1N1M1E1A1S0", "A1N1M1E1A1S1","A1N1M0E1A1S1")

#versions=("A1N0M0S0","A1N0M0S1","A1N0M0S2", "A1N7M0S0","A1N7M0S1", "A1N7M0S2") 
#versions=("A1N0M0S0","A1N0M0S1","A1N0M0S2", "A1N7M0S0",) 
#versions=("A1N0M0S0","A1N0M0S2", "A1N7M0S0","A1N7M0S2") 
#models=("A1N1")
#nametag="shuffling_NvsA"


## Fig. ACF_N_N4 a
#nametag="ACN"
#versions=("A2N1M0E1A1S1","A2N1M0E3A1S1","A2N1M0E1A0S1","A2N1M0E3A0S1")

#nametag="A2N1"
#versions=("A2N1M0E1A1S1",)



##nametag="ACN2"
##versions=("A2N1M0E1A1S1","A2N1M0E0A1S1","A2N1M0E1A0S1","A2N1M0E0A0S1")
##
##nametag="ACN3"
##versions=("A2N1M0E1A1S1","A2N1M0E2A1S1","A2N1M0E1A0S1","A2N1M0E2A0S1")

### Fig. ACF_N_N4 b
#nametag="ACN4" #E1A0 is unstable! now, all versions show an AC gap, including E3A0
#versions=("A1N1M0E1A1S1","A1N1M0E3A1S1","A1N1M0E3A0S1")


#nametag="A0"
##versions=("A1N1M1E1A1S0","A0N1M1E1A1S0",)
#versions=("A1N1M1E1A1S0","A0N1M1E1A1S0","A0N1M1E1A1S1", "A0N1M1E1A1S2")
#
#nametag="A3"
#versions=("A1N1M1E1A1S0","A3N1M1E1A1S0","A3N1M1E1A1S1", "A3N1M1E1A1S2")
#
#nametag="A0A3"
#versions=("A1N1M1E1A1S1","A3N1M1E1A1S1", "A0N1M1E1A1S1")

#nametag="A0N1M1"
#versions=("A1N1M1E1A1S1","A0N1M1E1A1S1",)

#nametag="A0N1M0"
#versions=("A1N1M0E1A1S1","A0N1M0E1A1S1",)

#nametag="A0Nx"
#versions=("A1N1M1E1A1S1","A0N1M1E1A1S1","A0N3M1E1A1S1")

#nametag="A0_M0"
#versions=("A1N1M0E1A1S1","A0N1M0E1A1S1","A0N3M0E1A1S1")

#nametag="A0N3"
##versions=("A1N1M0E1A1S1","A1N3M0E1A1S1","A0N3M0E1A1S1")
#versions=("A1N3M0E1A1S1","A0N3M0E1A1S1")

#nametag="A0N0"
#versions=("A0N3M0E1A1S1","A0N0M0E1A1S1")

#nametag="A1N0"
#versions=("A1N3M0E1A1S1","A1N0M0E1A1S1")


#############################

def plot_AC_all_fun(variable_name):
    plt.figure(dpi=dpi)
    plt.rcParams.update({'font.size': 30})
    fig, ax = plt.subplots(figsize=(10,10))
    variable=globals()['ACF_'+variable_name]
    #plt.title(variable_name)
    versionind=0

    for version in versions:
        #version=model+variant+'S'+str(shuffle)        
        model=version[0:4]
        shuffle=int(version[11])

        if versionind==0:
            color='#BB5566' #red
        elif versionind==1:
            color='#004488' # blue
        elif versionind==2:
            color='#DDAA33' #yellow
        elif versionind==3:
            color='black'

        if nametag=="ACN4":
            if versionind==0:
                color='#BB5566' #red
            elif versionind==1:
                color='#004488' # blue
            elif versionind==2:
                color='black'


        cellindind=0
        for cellind in cellindlist:
            if cellindind==0:
                style="-"
            elif cellindind==1:
                style="--"

            if nametag=="shuffling" or nametag=="AC_from_N_to_C" or \
            nametag=="AC_from_N_to_C2" or nametag=="AC_from_N_to_C5" or nametag=="shuffling_and_NoNdyn":
                if shuffle==0:
                    labelstr=model+' R0, P='+str(Pvec[cellind])+'mm/yr'
                elif shuffle==1:
                    labelstr=model+' R1, P='+str(Pvec[cellind])+'mm/yr'                            

#            elif nametag=="NoNdyn4":
#                if 
#                    labelstr=model+' M0, P='+str(Pvec[cellind])+'mm/yr'


            elif nametag[0:3]=="ACN":
                labelstr=version[6:10]+', P='+str(Pvec[cellind])+'mm/yr'
             
            else:
                labelstr=model+', P='+str(Pvec[cellind])+'mm/yr'

          #      labelstr=version+', P='+str(Pvec[cellind])+'mm/yr'
            
            
            plt.plot(variable[versionind,cellindind,:], linewidth=3, linestyle=style, color=color, label=labelstr)

            file_fig=exp + '_modelsandcells_'+nametag+'_ACF_'+variable_name+'.png'

            if nametag[0:3] == "PIK":
                ax.legend(fontsize=26)
            else:
                ax.legend(fontsize=18)
            
            plt.ylabel('autocorrelation')
            plt.xlabel('lag (years)')
            plt.ylim(-0.1,1)
            plt.xlim(0,lags)
            plt.savefig(file_fig, dpi=dpi)
            
            cellindind+=1
        versionind+=1


PFTind=PFT-1

Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]
Pvec=np.array(Pvec)                 # in mm/year


#Nmodels=len(models)
#Nvariants=len(variants)
#Nshuffles=len(shuffles)
#Nversions=Nmodels*Nvariants*Nshuffles
Nversions=len(versions)
ACF_Cperm2=np.zeros((Nversions,len(cellindlist),lags))
ACF_C=np.zeros((Nversions,len(cellindlist), lags))
ACF_N=np.zeros((Nversions,len(cellindlist), lags))

versionind=0
#for model in models:
#    for variant in variants:
#        for shuffle in shuffles:
for version in versions:

    #version=model+variant+'S'+str(shuffle)
    model=version[0:4]
    adj_switch=int(version[9])
    shuffle=int(version[11])

    YRS_spinup=2000
    
    if shuffle==0 or shuffle==3:
        YRS=10000
    else:
        YRS=52000        #  works thanks to shuffling


        #YRS=400       # TEST
        #YRS_spinup=100  # TEST

    
    #####
    YRS_LPJ=YRS
    if (YRS>10000):
        YRS_LPJ=10000


    
    ## drivers
    file = nc.Dataset(LPJdatapath+exp+'/pft_nind_beforemort.nc')
    N_beforemort = file['pft_nind']
    N_beforemort_LPJ=N_beforemort[0:YRS_LPJ,PFTind,index_PS,:]      
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_nind_beforeest.nc')
    N_beforeest = file['pft_nind']
    N_beforeest_LPJ=N_beforeest[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_nind_afterest.nc')
    N_afterest = file['pft_nind']
    N_afterest_LPJ=N_afterest[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_nind.nc')
    N_afteryear = file['nind']
    N_afteryear_LPJ=N_afteryear[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_bm_inc_carbon_beforeallo.nc',mode='r')
    bm = file.variables['bm_inc_carbon_beforeallo']
    bm_inc_LPJ=bm[0:YRS_LPJ,PFTind,index_PS,:]
    
    bm_inc_ind_LPJ=bm_inc_LPJ/N_beforemort_LPJ
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_wscal.nc')
    wscal_all = file['pft_wscal']
    wscal_LPJ=wscal_all[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_est.nc')
    est_read = file['establishment']
    est_LPJ=est_read[0:YRS_LPJ,PFTind,index_PS,:]

    file = nc.Dataset(LPJdatapath+exp+'/pft_mort.nc')
    mort_read = file['mortality']
    mort_LPJ=mort_read[0:YRS_LPJ,PFTind,index_PS,:]    
    
    ######## allocation fluxes
    # this is per ind, will be translated to perm2 using the fixed prescribed N
    file = nc.Dataset(LPJdatapath+exp+'/pft_tinc_ind_cleaf.nc')
    tinc_L = file['tinc_ind_cleaf']
    AL_LPJ=tinc_L[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_tinc_ind_croot.nc')
    tinc_R = file['tinc_ind_croot']
    AR_LPJ=tinc_R[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_tinc_ind_csapw.nc')
    tinc_S = file['tinc_ind_csapw']
    AS_LPJ=tinc_S[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_tinc_ind_chawo.nc')
    tinc_H = file['tinc_ind_chawo']
    AH_LPJ=tinc_H[0:YRS_LPJ,PFTind,index_PS,:]
    
    
    ## ind pools, used for version A2
    file = nc.Dataset(LPJdatapath+exp+'/pft_cleaf.nc')
    tinc_L = file['C_leaf_pft']
    L_LPJ=tinc_L[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_croot.nc')
    tinc_R = file['C_root_pft']
    R_LPJ=tinc_R[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_csapw.nc')
    tinc_S = file['C_sapwood_pft']
    S_LPJ=tinc_S[0:YRS_LPJ,PFTind,index_PS,:]
    
    file = nc.Dataset(LPJdatapath+exp+'/pft_chawo.nc')
    tinc_H = file['C_heartwood_pft']
    H_LPJ=tinc_H[0:YRS_LPJ,PFTind,index_PS,:]
        
    bm_inc_LPJ_mean=np.mean(bm_inc_LPJ, axis=0)    
    bm_inc_ind_LPJ_mean=np.mean(bm_inc_ind_LPJ, axis=0)
    wscal_LPJ_mean=np.mean(wscal_LPJ, axis=0)     
    N_afteryear_LPJ_mean=np.mean(N_afteryear_LPJ, axis=0)
    N_beforeest_LPJ_mean=np.mean(N_beforeest_LPJ, axis=0)
    N_afterest_LPJ_mean=np.mean(N_afterest_LPJ, axis=0)
    N_beforemort_LPJ_mean=np.mean(N_beforemort_LPJ, axis=0)
    
    AL_LPJ_mean=np.mean(AL_LPJ)
    AR_LPJ_mean=np.mean(AR_LPJ)        
    AS_LPJ_mean=np.mean(AS_LPJ)        
    AH_LPJ_mean=np.mean(AH_LPJ)
        
    A_LPJ=AL_LPJ+AR_LPJ+AS_LPJ+AH_LPJ
    A_LPJ_mean=np.mean(A_LPJ)
    
    L_LPJ_mean=np.mean(L_LPJ, axis=0)
    R_LPJ_mean=np.mean(R_LPJ, axis=0)
    S_LPJ_mean=np.mean(S_LPJ, axis=0)
    H_LPJ_mean=np.mean(H_LPJ, axis=0)
    
    ### timmean of allocated fractions in LPJ        
    #### note: A_LPJ is not identical to bm_inc_ind_LPJ
    # because bm_inc_ind is modified during allocation
    # also, there is the debt flux in addition to L, R, S, H
    # Here, for the fixed allocated fractions, use the sum of outputs
    # for L, R, S, H to get a sum = 1.
    
    ### prescribes allocated fluxes:    
    ## fractions are always fixed here
    # no temporal variability of fractions exists
    # but fracts still differ between locations (climate input)
    ALfrac_LPJ=np.mean(AL_LPJ/A_LPJ)
    ARfrac_LPJ=np.mean(AR_LPJ/A_LPJ)
    ASfrac_LPJ=np.mean(AS_LPJ/A_LPJ)
    AHfrac_LPJ=np.mean(AH_LPJ/A_LPJ)
            
    ### mean, std and AC of N dynamics
    Nmean_LPJ=np.mean(N_afteryear_LPJ, axis=0)
    Nstd_LPJ=np.std(N_afteryear_LPJ, axis=0)
    
    #NAC_LPJ=np.zeros(np.shape(N_afteryear_LPJ)[1])
    # NAC at lag 10 shall be = NAC1**10
#    for ind in range (0,np.shape(N_afteryear_LPJ)[1]):
#    #    NAC10_LPJ=ACF_fun(N_afteryear_LPJ[:,ind], 11)[10]
#    #    NAC_LPJ[ind]=np.power(NAC10_LPJ,1/10)
#        NAC30_LPJ=ACF_fun(N_afteryear_LPJ[:,ind], 31)[30]
#        NAC_LPJ[ind]=np.power(NAC30_LPJ,1/30)
    #sigma_noise=Nstd_LPJ*np.sqrt(1-NAC_LPJ)
    
    ## default parameters
    Nvarfrac=1
    varfrac=1
    mort_prescribed_param=-999
    est_prescribed=-999
    
    #if parameter=="fpc2":
    #    n_woody=2
    #else:
    #    n_woody=1
    n_woody=1
    
    par.set_fixed_parameters(PFT)
    fL=par.fL_default
    fR=par.fR_default
    fS=par.fS_default
    k_est=par.k_est_default     # same as est_max, maximum overall sapling establishment rate (indiv/m2)
    fpcmax=par.fpc_tree_max_default     # same as est_max, maximum overall sapling establishment rate (indiv/m2)
    
    
    ## set parameter values
    par.set_variable_parameters(Nparams, parameter, paramind)
    Nvarfrac=par.Nvarfrac
    varfrac=par.varfrac
    estfrac=par.estfrac
    AHfrac=par.AHfrac
    ASHfrac=par.ASHfrac
    mort_prescribed_param=par.mort_prescribed_param
    est_prescribed_param=par.est_prescribed_param
    fL=par.fL
    fR=par.fR
    fS=par.fS
    k_est=par.k_est
    n_woody=par.n_woody
    fpcmax=par.fpcmax
    param_name=par.param_name
    
    
    Nvarfrac_vec=par.Nvarfrac_vec
    varfrac_vec=par.varfrac_vec
    AHfrac_vec=par.AHfrac_vec
    ASHfrac_vec=par.ASHfrac_vec
    fL_vec=par.fL_vec
    fR_vec=par.fR_vec                    
    fS_vec=par.fS_vec
    mort_vec=par.mort_vec
    fpcmax_vec=par.fpcmax_vec
    n_woody_vec=par.n_woody_vec
    k_est_vec=par.k_est_vec
    est_vec=par.est_vec
    estfrac_vec=par.estfrac_vec
    
    
    if parameter=="AHfrac":
        ALfrac=ALfrac_LPJ
        ARfrac=ARfrac_LPJ
        ASfrac=1-ALfrac-ARfrac-AHfrac
    elif parameter=="ASHfrac":
        LofLR=ALfrac_LPJ/(ALfrac_LPJ+ARfrac_LPJ)
        SofSH=ASfrac_LPJ/(ASfrac_LPJ+AHfrac_LPJ)
        ALRfrac=1-ASHfrac
        ALfrac=ALRfrac*LofLR
        ARfrac=ALRfrac*(1-LofLR)
        ASfrac=ASHfrac*SofSH
        AHfrac=ASHfrac*(1-SofSH)
        
    else:
        ALfrac=ALfrac_LPJ
        ARfrac=ARfrac_LPJ
        ASfrac=ASfrac_LPJ
        AHfrac=AHfrac_LPJ
    
    if varfrac==1:
        bm_inc_input=bm_inc_LPJ
        bm_inc_ind_input=bm_inc_ind_LPJ
        wscal_input=wscal_LPJ
        
        AL_input=ALfrac*A_LPJ
        AR_input=ARfrac*A_LPJ
        AS_input=ASfrac*A_LPJ                
        AH_input=AHfrac*A_LPJ

        L_input=L_LPJ
        R_input=R_LPJ
        S_input=S_LPJ
        H_input=H_LPJ
        
    else:
        bm_inc_input=varfrac*(bm_inc_LPJ-bm_inc_LPJ_mean)+bm_inc_LPJ_mean
        bm_inc_ind_input=varfrac*(bm_inc_ind_LPJ-bm_inc_ind_LPJ_mean)+bm_inc_ind_LPJ_mean
        wscal_input=varfrac*(wscal_LPJ-wscal_LPJ_mean)+wscal_LPJ_mean        

        AL_input=ALfrac*(varfrac*(A_LPJ-A_LPJ_mean)+A_LPJ_mean)
        AR_input=ARfrac*(varfrac*(A_LPJ-A_LPJ_mean)+A_LPJ_mean)       
        AS_input=ASfrac*(varfrac*(A_LPJ-A_LPJ_mean)+A_LPJ_mean)
        AH_input=AHfrac*(varfrac*(A_LPJ-A_LPJ_mean)+A_LPJ_mean)

        L_input=varfrac*(L_LPJ-L_LPJ_mean)+L_LPJ_mean
        R_input=varfrac*(R_LPJ-R_LPJ_mean)+R_LPJ_mean       
        S_input=varfrac*(S_LPJ-S_LPJ_mean)+S_LPJ_mean
        H_input=varfrac*(H_LPJ-H_LPJ_mean)+H_LPJ_mean

    
    if Nvarfrac==1:
        N_afteryear_input=N_afteryear_LPJ
        N_beforeest_input=N_beforeest_LPJ    
        N_afterest_input=N_afterest_LPJ
        N_beforemort_input=N_beforemort_LPJ    
    
    else:
        N_afteryear_input=Nvarfrac*(N_afteryear_LPJ-N_afteryear_LPJ_mean)+N_afteryear_LPJ_mean    
        N_beforeest_input=Nvarfrac*(N_beforeest_LPJ-N_beforeest_LPJ_mean)+N_beforeest_LPJ_mean    
        N_afterest_input=Nvarfrac*(N_afterest_LPJ-N_afterest_LPJ_mean)+N_afterest_LPJ_mean    
        N_beforemort_input=Nvarfrac*(N_beforemort_LPJ-N_beforemort_LPJ_mean)+N_beforemort_LPJ_mean    
    
    


    cellindind=0
    for cellind in cellindlist:

        
        ####### initialise arrays
        L_yearly=np.zeros(YRS)
        R_yearly=np.zeros(YRS)
        S_yearly=np.zeros(YRS)
        H_yearly=np.zeros(YRS)
        D_yearly=np.zeros(YRS)
        N_yearly=np.zeros(YRS)

        fpc_yearly=np.zeros(YRS)
        C_yearly=np.zeros(YRS)
        Cperm2_yearly=np.zeros(YRS)
        
#        AL_yearly=np.zeros(YRS)
#        AR_yearly=np.zeros(YRS)
#        AS_yearly=np.zeros(YRS)
#        AH_yearly=np.zeros(YRS)
      
    
        # inicond in case last cell ended up in a mess
        L, R, S, H, D, N, fpc, height = 9057.575, 9157.507, 107325.516, 140601.61, 0.0, 0.06495886, 0.95, 17.915508
    
           
        if (version[4:6]=="M0" and not parameter=="mort"):
            mort_prescribed=0.015                    
        elif (version[4:6]=="M2" and not parameter=="mort"):
            mort_prescribed=np.mean(mort_LPJ[:,cellind])
        else:
            mort_prescribed=par.mort_prescribed_param

        if (version[2:4]=="N8" and not parameter=="est"):
            est_prescribed=np.mean(est_LPJ[:,cellind])
        else:
            est_prescribed=par.est_prescribed_param

        

        for year in range(1,YRS):
    

            if shuffle==0:   # no shuffling in any input
                year_A=year
                year_N=year
                year_M=year
                year_E=year
    
            elif shuffle==1:  # shuffle all inputs the same way (destroys AC, but keeps corr)
                year_A=random.randint(1, YRS_LPJ-1)
                year_N=year_A
                year_M=year_A
                year_E=year_A
                
            elif shuffle==2:  # shuffle all inputs differently (destroy corr)
                year_A=random.randint(1, YRS_LPJ-1)
                year_N=random.randint(1, YRS_LPJ-1)
                year_M=random.randint(1, YRS_LPJ-1)
                year_E=random.randint(1, YRS_LPJ-1)
                
            elif shuffle==3:  # shuffle allocated C, but keep N ordered
                year_A=random.randint(1, YRS_LPJ-1)
                year_N=year
                year_M=year
                year_E=year
            
                                                   
            
            if (version[0:2]=="A0"):
                ## use fixed fractions and distribute across these fractions: 
                ALperm2=AL_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                ARperm2=AR_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                ASperm2=AS_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                AHperm2=AH_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                
            elif (version[0:2]=="A3"):
                ## randomised allocated fractions: 
                ## read allocated fluxes from LPJ (fractions can vary)
                 # just as a sanity check: S0 (no shuffle) gives identical results to version A1
                ALperm2=AL_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                ARperm2=AR_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                ASperm2=AS_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                AHperm2=AH_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]


            if (version[0:2]=="A2"):
                L=L_input[year_A, cellind]
                R=R_input[year_A, cellind]
                S=S_input[year_A, cellind]
                H=H_input[year_A, cellind]
                height=0
    
            if (version[2:4]=="N0"):
                N=N_afteryear_input[year_N, cellind]
                
            elif (version[2:4]=="N2"):
                Nold=N_beforeest_input[year_N, cellind]
                Nnew=N_afterest_input[year_N, cellind]
                est=Nnew-Nold
                N=N_afteryear_input[year_N, cellind]
                estvsN=est/Nold
            elif (version[2:4]=="N3"):
                N=Nmean_LPJ[cellind]
            #elif (version[2:4]=="N7"):   # AR1 process
            #    if shuffle==0:
            #        N = NAC_LPJ[cellind]*(N-Nmean_LPJ[cellind]) + random.gauss(Nmean_LPJ[cellind],sigma_noise[cellind])
            #    else:
            #        N = random.gauss(Nmean_LPJ[cellind],Nstd_LPJ[cellind])
    
    
            if parameter=="mort": 
                mort_prescribed=par.mort_prescribed_param
            elif (version[4:6]=="M0"):
                mort_prescribed=0.015
            elif version[4:6]=="M1":
                mort_prescribed=-999
            elif (version[4:6]=="M2"):
                mort_prescribed=np.mean(mort_LPJ[:,cellind])
            elif (version[4:6]=="M3"):
                mort_prescribed=mort_LPJ[year_M,cellind]                
            
            if parameter=="est":
                est_prescribed=par.est_prescribed_param
            elif (version[6:8]=="E0"):
                est_prescribed=0.01
            elif version[6:8]=="E1":
                est_prescribed=-999
            elif (version[6:8]=="E2"):
                est_prescribed=np.mean(est_LPJ[:,cellind])
            elif (version[6:8]=="E3"):
                est_prescribed=est_LPJ[year_E,cellind]                


            ### run the model
            ## two options, which are theoretically identical if input is not shuffled, but which matter for numerical stability:
            bm_forcing = bm_inc_ind_input[year_A, cellind]  # this is the default mode; passes bm_inc per individual from LPJ
            #bm_forcing = bm_inc_input[year_A, cellind]/N  # alternative mode; passes bm_inc per m2 from LPJ
    
            if model=='A1N1':
                L, R, S, H, D, N, fpc, height, AL, AR, AS, AH, drought1, drought2, \
                    adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest = \
                       LPJ_functions.LPJ_Cbalance_A1N1(L, R, S, H, \
                          D, N, height, bm_inc_ind_input[year_A, cellind], \
                          wscal_input[year_A, cellind], fL, fR, fS, Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch)
                ## not used in paper (E always shuffled or interactive there):
                ## in case E and M are fixed, there is no variability in N. In this case, add noise:
                if (parameter=="est" or version[6:8]=="E0" or version[6:8]=="E2") and (parameter=="mort" or version[4:6]=="M0" or version[4:6]=="M2"):
                    N=N+random.gauss(0,Nstd_LPJ[cellind])
                
                    
            elif model=='A0N1' or model=='A3N1':
                L, R, S, H, N, fpc, height, AL, AR, AS, AH, \
                    adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest = \
                    LPJ_functions.LPJ_Cbalance_A0N1(L, R, S, H, \
                       N, height, ALperm2, ARperm2, ASperm2, AHperm2, fL, fR, fS, Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch)
    
            elif model=='A1N2':
                L, R, S, H, D, height, AL, AR, AS, AH, drought1, drought2 = \
                    LPJ_functions.LPJ_Cbalance_A1N2(L, R, S, H, \
                       D, Nold, est, height, bm_inc_ind_input[year_A, cellind], \
                       wscal_input[year_A, cellind], fL, fR, fS)
    
            elif model=='A1N0' or model=='A1N7' or model=='A1N3':  # prescribed N, use from LPJ, AR1 or constant
                L, R, S, H, D, height, AL, AR, AS, AH, drought1, drought2, mort, bm_delta, turnover = \
                    LPJ_functions.LPJ_Cbalance_A1N0(L, R, S, H, \
                       D, height, bm_inc_ind_input[year_A, cellind], \
                       wscal_input[year_A, cellind], fL, fR, fS, mort_prescribed)


### here: feed in A perm2
#                elif model=='A0N0' or model=='A0N3' or model=='A0N7':
#                    L, R, S, H, mort, bm_delta, turnover = LPJ_functions.LPJ_Cbalance_A0N0(L, R, S, H, N, \
#                      ALperm2, ARperm2, ASperm2, AHperm2, fL, fR, fS, mort_prescribed)
#
#                    ## only for diagnostic purpose:                            
#                    AL=ALperm2/N
#                    AR=ARperm2/N
#                    AS=ASperm2/N
#                    AH=AHperm2/N
    
## feed in A per ind
            elif model=='A0N0' or model=='A0N3' or model=='A0N7':
                
                ## only for diagnostic purpose:                            
                AL=AL_input[year_A, cellind]
                AR=AR_input[year_A, cellind]
                AS=AS_input[year_A, cellind]
                AH=AH_input[year_A, cellind]
                
                L, R, S, H, mort, bm_delta, turnover = LPJ_functions.LPJ_Cbalance_A0N0(L, R, S, H, \
                  AL, AR, AS, AH, fL, fR, fS, mort_prescribed)

    
    
            elif model=='A2N1':                         
                N, fpc, adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest = \
                   LPJ_functions.LPJ_Cbalance_A2N1(L, R, S, H, D, AL, AR, AS, AH, N, fL, fR, fS, \
                      Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch)
        
                ## in case E and M are fixed, there is no variability in N. In this case, add noise:
                if (parameter=="est" or version[6:8]=="E0" or version[6:8]=="E2") and (parameter=="mort" or version[4:6]=="M0" or version[4:6]=="M2"):
                    N=N+random.gauss(0,Nstd_LPJ[cellind])        

            else:
                print("model version not found")
        
            # save this year's output
            L_yearly[year]=L
            R_yearly[year]=R
            S_yearly[year]=S
            H_yearly[year]=H
            Cperm2_yearly[year]=(L+R+S+H)*N
            C_yearly[year]=L+R+S+H
            N_yearly[year]=N
    

        ACF_Cperm2[versionind,cellindind,:]=ACF_fun(Cperm2_yearly[YRS_spinup:], lags)
        ACF_C[versionind,cellindind,:]=ACF_fun(C_yearly[YRS_spinup:], lags)
        ACF_N[versionind,cellindind,:]=ACF_fun(N_yearly[YRS_spinup:], lags)

            
        if plot_AC_single==1:
            plt.figure(dpi=dpi)
            plt.rcParams.update({'font.size': 30})
            fig, ax = plt.subplots(figsize=(10,10))
            plt.plot(ACF_Cperm2[versionind,cellindind,:], color="black", linewidth=3, linestyle="-")
            file_fig=exp + '_' + version +'_ACF_Cperm2_cellind'+str(cellind)+'.png'
            plt.ylabel('ACF')
            plt.xlabel('lag (years)')
            plt.savefig(file_fig, dpi=dpi)
            plt.xlim(0,lags)
       
        cellindind=cellindind+1    
    


    if plot_AC_cellpairs==1:
        plt.figure(dpi=dpi)
        plt.rcParams.update({'font.size': 30})
        fig, ax = plt.subplots(figsize=(10,10))
        plt.plot(ACF_Cperm2[versionind,0,:], color="black", linewidth=3, linestyle="-", label="wet")
        plt.plot(ACF_Cperm2[versionind,1,:], color="black", linewidth=3, linestyle="--", label="dry")
        #file_fig="model"+model+"_case"+str(case)+"_AFC_Cperm2.png"
        file_fig=exp + '_' + version +'_ACF_Cperm2_cellpair.png'
        ax.legend()
        plt.ylabel('ACF')
        plt.xlabel('lag (years)')
        plt.ylim(0,1)
        plt.xlim(0,lags)

        plt.savefig(file_fig, dpi=dpi)

    versionind=versionind+1


if plot_AC_versions==1:
    cellindind=0
    for cellind in cellindlist:
        plt.figure(dpi=dpi)
        plt.rcParams.update({'font.size': 30})
        fig, ax = plt.subplots(figsize=(10,10))
        versionind=0
        style='-'            
        for version in versions:

            #version=model+variant+'S'+str(shuffle)        
            plt.plot(ACF_Cperm2[versionind,cellindind,:], linewidth=3, linestyle=style, label=version)
            versionind+=1
        file_fig=exp + '_modelpair_'+parameter+'_ACF_Cperm2_cellind'+str(cellind)+'.png'
        ax.legend()
        plt.ylabel('ACF')
        plt.xlabel('lag (years)')
        plt.ylim(0,1)
        plt.xlim(0,lags)
        plt.savefig(file_fig, dpi=dpi)
        cellindind+=1


if plot_AC_all==1:
    if nametag[0:3]=="ACN":
        plot_AC_all_fun("N")
        plot_AC_all_fun("Cperm2")  ## comment out
    elif nametag[0:2]=="A0" or nametag=="A0_M0":
        plot_AC_all_fun("Cperm2")
    else:
        plot_AC_all_fun("Cperm2")
        #plot_AC_all_fun("C")
        #plot_AC_all_fun("N")
