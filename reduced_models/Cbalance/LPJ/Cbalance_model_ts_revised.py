#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 17:05:13 2023

@author: bathiany
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
sys.path.insert(1, '/home/bathiany/Projects/Vegetation_resilience_indicators/reduced_models/Cbalance/functions')
import LPJ_functions
import par
#from ACF import ACF_fun
from sklearn.linear_model import LinearRegression
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
np.seterr(invalid='ignore')

LPJdatapath='/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/'

########## parameters
dpi=500
exp="sba0084"      # exp from which data is read in
YRS_spinup=2000
PFT=1

cellindlist=[2,]
# wet cell, 2000 mm, has cellind 10,
# dry cell, 400 mm, has cellind 2

index_PS=0

#models=("A1N1",)
### effect of M, E, A on AC of N (A2)
#variants=("M1E1A1","M0E0A1","M2E2A1","M3E3A1",\
#          "M0E1A1","M1E0A1",\
#          "M2E0A1","M0E2A1",\
#          "M2E1A1","M1E2A1",\
#          "M3E1A1","M1E3A1",\
#          "M0E3A1","M3E0A1",\
#          "M2E3A1","M3E2A1",\
#          "M1E1A0","M0E0A0","M2E2A0","M3E3A0",\
#          "M0E1A0","M1E0A0",\
#          "M2E0A0","M0E2A0",\
#          "M2E1A0","M1E2A0",\
#          "M3E1A0","M1E3A0",\
#          "M0E3A0","M3E0A0",\
#          "M2E3A0","M3E2A0")


#versions=("A1N1M1E1A1S0",) # default, unshuffled for eval
#versions=("A1N1M1E1A1S1",) # default, shuffled

#versions=("A3N1M1E1A1S1",) # shuffled allocated fluxes 
#versions=("A0N1M1E1A1S1",)  # fixed allocated fractions

#versions=("A5N1M1E1A1S1",)  # like default but forced with ind-based fluxes

## when plotting A fluxes and their fractions based on ind-based fluxes,
# this distorts the scatter cloud; e.g. Cass becomes larger!
# the distrib is only conserved when using the perm2-based values
# because these are the ones used as input to drive the model

# could this be problematic when modifying the allocation?
# if model has larger N than LPJ sba0044, the same Aperm2 is distributed over more trees
# => each one gets less C... => in this sense, C is not conserved for each tree!
# This might actually generate oscillations => large AC??


## strange, fails
#variants=("M1E1A0", "M0E1A0", "M2E1A0", "M3E1A0")

# scatter, stability
#versions=("A1N1M0E1A1S1","A2N1M0E1A1S1",)

#versions=("A1N1M1E1A1S1",)  
## E+M+A has positive slope in the default model, because of mortality
## single line is not a good fit, the tendencies are bimodal.

# allocation scatter plots
versions=("A1N1M0E1A1S1",)


#models=("A1N1",)
#variants=("M1E1A1", "M1E1A0", "M1E3A1", "M1E3A0")

#versions=("A1N1M1E1A1S1",)

#models=("A1N1",)
### effect of M, E, A on AC of N (A2)
#variants=("M1E3A0","M1E1A0","M1E3A1","M1E2A0","M1E2A1")

## A0
#versions=("A1N1M1E1A1S0","A0N1M1E1A1S0",)

#versions=("A1N1M1E1A1S0",) 
#versions=("A0N1M1E1A1S0",)

#versions=("A0N1M1E1A1S1",)
#versions=("A0N1M1E1A1S1",)
#versions=("A0N0M1E1A1S1",)
#versions=("A1N0M1E1A1S1",)

## allocation & variance
#versions=("A1N1M0E1A1S1",)

### use varfrac with A0:
#versions=("A0N3M1E1A1S1",)



## what role do M, E, adj play; in Ndyn only and in coupled version

# default model
#parameter="varfrac"
#paramind=5 # default
#paramind=6 #1.5
#paramind=2 #0.1

parameter="Nvarfrac"
paramind=5 #default

## ts visual inspection
#YRS_spinup=2000
#YRS=1000+YRS_spinup

### # N and allocation scatter plots
YRS_spinup=2000
YRS=4000  

### eval ts
#YRS_spinup=2000
#YRS=250+YRS_spinup

## recovery
#YRS_spinup=1
#YRS=100


Nini=1


######
Nparams=7
#par.set_fixed_parameters(PFT)


YRS_LPJ=YRS
if (YRS>10000):
    YRS_LPJ=10000

def fig_ts(ts, ts_name, paramind, cellind, version, shuffle, parameter):
    plt.figure(dpi=dpi)
    plt.rcParams.update({'font.size': 22})
    fig, ax = plt.subplots(figsize=(8,6))
    plt.plot(ts, color="black", linewidth=1.5)

    if ts_name=="Cperm2":
        plt.ylabel("$C \; [g/m^2]$")
    #    filename_AC=exp + '_' + version + 'S'+str(shuffle)+'_'+parameter+'_vec_AC10_Cperm2_diag.npy'
        #data=np.load(filename_AC)
                                            
    elif ts_name=="C":
        plt.ylabel("$C \; [g/ind]$")
    elif ts_name=="N" or ts_name=="N_LPJ":
        plt.ylabel("$N \; [ind/m^2]$")
    elif ts_name=="fpc":
        plt.ylabel("$fpc \; [m^2/m^2]$")
    elif ts_name=="est":
        plt.ylabel("$establishment \; [ind/m^2yr]$")
    elif ts_name=="mort":
        plt.ylabel("$mortality rate \; [1/yr]$")
    else:
        plt.ylabel(ts_name)

    plt.xlabel('$year$')
    plt.title(version+', paramind '+ str(paramind)+', cellind '+str(cellind), fontsize=18)
    file_fig=exp+'_'+version+'_'+parameter+'_paramind'+ str(paramind)+'_cellind'+str(cellind)+'_'+ts_name+'_ts.png'
    plt.savefig(file_fig)     
    

def fig_ts_eval(ts1, ts2, ts_name, paramind, cellind, version, shuffle, parameter):
    plt.figure(dpi=dpi)
    plt.rcParams.update({'font.size': 22})
    fig, ax = plt.subplots(figsize=(8,6))

    plt.plot(ts1, color="black", linewidth=2, label="reduced model")
    plt.plot(ts2, color="red", linewidth=0.5, label="LPJ")

    #ax.legend(loc='upper left')


    if ts_name=="Cperm2":
        plt.ylabel(r'$\tilde{C} \; [g/m^2]$')                                            
    elif ts_name=="C":
        plt.ylabel("$C \; [g/ind]$")
    elif ts_name=="N" or ts_name=="N_LPJ":
        plt.ylabel("$N \; [ind/m^2]$")
    elif ts_name=="fpc":
        plt.ylabel("$fpc \; [m^2/m^2]$")
    elif ts_name=="est":
        plt.ylabel("$establishment \; [ind/m^2yr]$")
    elif ts_name=="mort":
        plt.ylabel("$mortality rate \; [1/yr]$")
    else:
        plt.ylabel(ts_name)

    plt.xlabel('$year$')
    plt.title(ts_name + ', ' + version+', paramind '+ str(paramind)+', cellind '+str(cellind), fontsize=18)
    file_fig=exp+'_'+version+'_'+parameter+'_paramind'+ str(paramind)+'_cellind'+str(cellind)+'_'+ts_name+'_ts_eval.png'
    plt.savefig(file_fig)     
    

def fig_N_scatter(ts1,ts2,ts1_name,ts2_name):
    plt.figure(dpi=dpi)
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(8,6))
    
    plt.scatter(ts1, ts2, s=6)
#    plt.xlabel(ts1_name)
#    plt.ylabel(ts2_name)
    plt.xlabel(ts1_name + ' [ind/m²]')
    plt.ylabel(ts2_name + ' [ind/m²yr]')

    regressor = LinearRegression() 
    regressor.fit(ts1.reshape((-1, 1)), ts2)

    y_pred = regressor.predict(ts1.reshape((-1, 1)))
    plt.plot(ts1, y_pred,color='r')

    if ts2_name=="E+M+A":
        plt.axhline(y = 0, color = 'k', linestyle='--', linewidth=1)   # dN/dt = 0
        plt.axvline(x = np.mean(ts1), color = 'k', linestyle='--', linewidth=1) # N mean

    plt.title(version+', cellind '+str(cellind), fontsize=18)
    file_fig=exp+'_'+version+'_'+parameter+'_paramind'+ str(paramind)+'_cellind'+str(cellind)+'_scatter_'+ts1_name+'_'+ts2_name+'.png'
    plt.xticks(rotation=90)
    plt.savefig(file_fig)
    plt.show()



def fig_A_scatter(var1, var2, var3, var1_name, var2_name, var3_name):    
    plt.figure(dpi=dpi)
    
    plt.rcParams.update({'font.size': 24})
    fig, ax = plt.subplots(figsize=(8,6))
    #var1=globals()[var1_name + '_yearly']
    #var2=globals()[var2_name + '_yearly']
    #var3=globals()[var3_name + '_yearly']
    plt.scatter(var1, var2, s=6, c=var3)   #s=6
    plt.xlabel(var1_name)
    plt.ylabel(var2_name)
    if var1_name=='A':
        plt.xlabel('$C_{ass} \: [gC/yr]$')

    if var2_name=='AS':
        plt.ylabel('$A_S \: [gC/yr]$')
    elif var2_name=='ALR':
        plt.ylabel('$A_L + A_R \: [gC/yr]$')        
    elif var2_name=='ASH':
        plt.ylabel('$A_S + A_H \: [gC/yr]$')        
    elif var2_name=='ASHfrac':
        plt.ylabel('$(A_S + A_H)/C_{ass}$')   
        #plt.ylim(0.2,0.5)

    #plt.title()

    ## time mean:
    var1mean=np.mean(var1)
    var2mean=np.mean(var2)
    plt.axvline(x = var1mean, color = 'r', linestyle='--', linewidth=2)    
    plt.axhline(y = var2mean, color = 'r', linestyle='--', linewidth=2)

    #file_fig='allocation_'+var1_name+ '_vs_'+var2_name+'_scatter_gridcell'+str(gridcell)+'_varfrac'+str(varfrac)+'.png'
    file_fig=exp+'_'+version+'_'+parameter+'_paramind'+ str(paramind)+'_cellind'+str(cellind)+'_allocation_scatter_'+var1_name+'_'+var2_name+'.png'

    plt.savefig(file_fig)
    print(file_fig)
    plt.show()



Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]
Pvec=np.array(Pvec)                 # in mm/year

PFTind=PFT-1


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


file = nc.Dataset(LPJdatapath+exp+'/fpc.nc')
fpc_read = file['FPC']
fpc_LPJ=fpc_read[0:YRS_LPJ,PFTind,index_PS,:]   

    
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
    

file = nc.Dataset(LPJdatapath+exp+'/pft_vegc.nc')
vegc_H = file['vegC_pft']
Cperm2_LPJ=vegc_H[0:YRS_LPJ,PFTind,index_PS,:]



#print(Cperm2_LPJ[YRS_spinup:,10])
#print(Cperm2_LPJ[:,10])

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
## NAC at lag 10 shall be = NAC1**10
#for ind in range (0,np.shape(N_afteryear_LPJ)[1]):
##    NAC10_LPJ=ACF_fun(N_afteryear_LPJ[:,ind], 11)[10]
##    NAC_LPJ[ind]=np.power(NAC10_LPJ,1/10)
#    NAC30_LPJ=ACF_fun(N_afteryear_LPJ[:,ind], 31)[30]
#    NAC_LPJ[ind]=np.power(NAC30_LPJ,1/30)
#sigma_noise=Nstd_LPJ*np.sqrt(1-NAC_LPJ)

## default parameters
Nvarfrac=1
varfrac=1


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
#mort_prescribed_param=par.mort_prescribed_param
#est_prescribed_param=par.est_prescribed_param
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
    # affects how much transfer happens from S to H, rest is fixed.
    # ALfrac and ARfrac are unchanged, but more AH reduces AS.                                        
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
    

#for model in models:
#    for variant in variants:
#        version=model+variant+'S'+str(shuffle)
#        adj_switch=int(variant[5])
for version in versions:

        model=version[0:4]
        adj_switch=int(version[9])
        shuffle=int(version[11])
        #############################    
    
        
        for cellind in cellindlist:
                        
            # diagnosis           
            mort_yearly=np.zeros(YRS)
            Nmort_yearly=np.zeros(YRS)

            bm_delta_yearly=np.zeros(YRS)
            turnover_yearly=np.zeros(YRS)
            est_yearly=np.zeros(YRS)
            adj_yearly=np.zeros(YRS)
            adjN_yearly=np.zeros(YRS)
            adjNrel_yearly=np.zeros(YRS)
            drought1_yearly=np.zeros(YRS)
            drought2_yearly=np.zeros(YRS)
            drought_yearly=np.zeros(YRS)

            C_yearly=np.zeros(YRS)
            Cperm2_yearly=np.zeros(YRS)
            LSperm2_yearly=np.zeros(YRS)
        
            height_yearly=np.zeros(YRS)
        
            fpc_beforeest_yearly=np.zeros(YRS)
            N_beforeest_yearly=np.zeros(YRS)
        
            estvsN_yearly=np.zeros(YRS)
        
            N_yearly=np.zeros(YRS)
            fpc_yearly=np.zeros(YRS)
           
            L_yearly=np.zeros(YRS)
            R_yearly=np.zeros(YRS)
            S_yearly=np.zeros(YRS)
            H_yearly=np.zeros(YRS)
        
            AL_yearly=np.zeros(YRS)
            AR_yearly=np.zeros(YRS)
            AS_yearly=np.zeros(YRS)
            AH_yearly=np.zeros(YRS)
            A_yearly=np.zeros(YRS)

            ALperm2_yearly=np.zeros(YRS)
            ARperm2_yearly=np.zeros(YRS)
            ASperm2_yearly=np.zeros(YRS)
            AHperm2_yearly=np.zeros(YRS)
            Aperm2_yearly=np.zeros(YRS)

            SHfrac_yearly=np.zeros(YRS)
            ASHfrac_yearly=np.zeros(YRS)
            Hfrac_yearly=np.zeros(YRS)
            AHfrac_yearly=np.zeros(YRS)
        
            ## for eval      
            AL_LPJ_yearly=np.zeros(YRS)
                
            Czero_yearly=np.zeros(YRS)
                
            # inicond in case last cell ended up in a mess
            L, R, S, H, D, N, fpc, height = 9057.575, 9157.507, 107325.516, 140601.61, 0.0, 0.06495886, 0.95, 17.915508
        
            #N=0 # inicond
            L, R, S, H, D, N, fpc, height = 0, 0, 0, 0, 0, 0, 0, 0
            N=Nini

            # initialise from LPJ...
            ######## inicond
            L_yearly[0]=L_LPJ[0, cellind]
            R_yearly[0]=R_LPJ[0, cellind]
            S_yearly[0]=S_LPJ[0, cellind]
            H_yearly[0]=H_LPJ[0, cellind]
            N_yearly[0]=N_afteryear_LPJ[0, cellind]
            fpc_yearly[0]=fpc_LPJ[0, cellind]
            #### initialise state variables
            L=L_yearly[0]
            R=R_yearly[0]
            S=S_yearly[0]
            H=H_yearly[0]
            N=N_yearly[0]
            fpc=fpc_yearly[0]


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
                
                
                ## set outputs to 0, then models will overwrite those that they compute
                drought1=0
                drought2=0
                adjN=0
                adjNrel=0
                adj_count=0
                mort=0
                bm_delta=0
                turnover=0                    
                est=0
                estvsN=0
                adj_count=0
                fpc_beforeest=0
                N_beforeest=0

                if (version[0:2]=="A0"):
                    ## use fixed fractions and distribute across these fractions: 
                    ALperm2=AL_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    ARperm2=AR_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    ASperm2=AS_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    AHperm2=AH_input[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    Aperm2=ALperm2+ARperm2+ASperm2+AHperm2

                elif (version[0:2]=="A3"):
                    ## randomised allocated fractions: 
                    ## read allocated fluxes from LPJ (fractions can vary)
                     # just as a sanity check: S0 (no shuffle) gives identical results to version A1
                    ALperm2=AL_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    ARperm2=AR_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    ASperm2=AS_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    AHperm2=AH_LPJ[year_A, cellind]*N_afteryear_LPJ[year_A-1, cellind]
                    Aperm2=ALperm2+ARperm2+ASperm2+AHperm2

                elif (version[0:2]=="A4"):    ## like A0, but forced with per indiv fluxes - unstable, develops oscillations           
                    AL=AL_LPJ[year_A, cellind]
                    AR=AR_LPJ[year_A, cellind]
                    AS=AS_LPJ[year_A, cellind]
                    AH=AH_LPJ[year_A, cellind]
                    A=AL+AR+AS+AH

                elif (version[0:2]=="A5"):  # like A3 but with ind-based forcing
                    AL=AL_input[year_A, cellind]
                    AR=AR_input[year_A, cellind]
                    AS=AS_input[year_A, cellind]
                    AH=AH_input[year_A, cellind]
                    A=AL+AR+AS+AH

                if version[0:2]=="A2":
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
#                elif (version[2:4]=="N7"):   # AR1 process
#                    if shuffle==0:
#                        N = NAC_LPJ[cellind]*(N-Nmean_LPJ[cellind]) + random.gauss(Nmean_LPJ[cellind],sigma_noise[cellind])
#                    else:
#                        N = random.gauss(Nmean_LPJ[cellind],Nstd_LPJ[cellind])
        
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

                ## two options, which are theoretically identical if input is not shuffled, but which matter for numerical stability:
                bm_forcing = bm_inc_ind_input[year_A, cellind]  # this is the default mode; passes bm_inc per individual from LPJ
                #bm_forcing = bm_inc_input[year_A, cellind]/N  # alternative mode; passes bm_inc per m2 from LPJ

                ### run the model
                if model=='A1N1':
                                        
                    L, R, S, H, D, N, fpc, height, AL, AR, AS, AH, drought1, drought2, \
                        adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest = \
                           LPJ_functions.LPJ_Cbalance_A1N1(L, R, S, H, \
                              D, N, height, bm_forcing, \
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
        
                elif model=='A4N1' or model=='A5N1':   
                    L, R, S, H, N, fpc, height, ALperm2, ARperm2, ASperm2, AHperm2, \
                        adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest = \
                        LPJ_functions.LPJ_Cbalance_A4N1(L, R, S, H, \
                           N, height, AL, AR, AS, AH, fL, fR, fS, Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch)
                        
                elif model=='A1N2':
                    L, R, S, H, D, height, AL, AR, AS, AH, drought1, drought2 = \
                        LPJ_functions.LPJ_Cbalance_A1N2(L, R, S, H, \
                           D, Nold, est, height, bm_forcing, \
                           wscal_input[year_A, cellind], fL, fR, fS)
        
                elif model=='A1N0' or model=='A1N7' or model=='A1N3':  # prescribed N, use from LPJ, AR1 or constant
                    L, R, S, H, D, height, AL, AR, AS, AH, drought1, drought2, mort, bm_delta, turnover = \
                        LPJ_functions.LPJ_Cbalance_A1N0(L, R, S, H, \
                           D, height, bm_forcing, \
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
                LSperm2_yearly[year]=(L+S)*N
                C_yearly[year]=L+R+S+H
                AL_yearly[year]=AL
                AR_yearly[year]=AR
                AS_yearly[year]=AS
                AH_yearly[year]=AH           
                A_yearly[year]=AL+AR+AS+AH
        
                if version[0:2]=="A0" or version[0:2]=="A3" or version[0:2]=="A4":
                    ALperm2_yearly[year]=ALperm2
                    ARperm2_yearly[year]=ARperm2
                    ASperm2_yearly[year]=ASperm2
                    AHperm2_yearly[year]=AHperm2      
                    Aperm2_yearly[year]=Aperm2

                N_yearly[year]=N
                fpc_yearly[year]=fpc
        
                fpc_beforeest_yearly[year]=fpc_beforeest
                N_beforeest_yearly[year]=N_beforeest
        
                height_yearly[year]=height
        
                drought1_yearly[year]=drought1
                drought2_yearly[year]=drought2
                drought_yearly[year]=drought1+drought2

                adj_yearly[year]=adj_count
                adjN_yearly[year]=adjN
                adjNrel_yearly[year]=adjNrel
                mort_yearly[year]=mort
                Nmort_yearly[year]=Nmort
                bm_delta_yearly[year]=bm_delta
                turnover_yearly[year]=turnover
                est_yearly[year]=est
                estvsN_yearly[year]=estvsN

                SHfrac_yearly[year]=(S+H)/(L+R+S+H)
                ASHfrac_yearly[year]=(AS+AH)/(AL+AR+AS+AH)
                Hfrac_yearly[year]=H/(L+R+S+H)
                AHfrac_yearly[year]=AH/(AL+AR+AS+AH)
        
                if L<1e-9 or R<1e-9 or S<1e-9 or H<1e-9:
                    Czero_yearly[year]=1
        
                ## for comparison: LPJ, for example to compare the 
                # allocated fluxes to the ones allocated in LPJ for the same
                # NPP and soil moisture
                # works as eval of allo module
                #AL_LPJ_yearly[year]=AL_LPJ[year_A, cellind]
        
        
                #### end of time loop
                #Cperm2_anomaly_bw40, Cperm2_smoothed_bw40 = ts_separate(Cperm2_yearly[YRS_spinup:], bw=40)
                #Cperm2_anomaly_bw100, Cperm2_smoothed_bw100 = ts_separate(Cperm2_yearly[YRS_spinup:], bw=100)
 
 
           # time series for sanity checks by eye
            #if asym_diag[paramind,cellind]>4 or ACFmin_diag[paramind,cellind]<-0.6:
            #if ( cellind==cellind1_plot or cellind==cellind2_plot) and paramind==paramind_plot:
#            fig_ts(Cperm2_yearly[YRS_spinup:], 'Cperm2', paramind, cellind, version, shuffle, parameter)
            #fig_ts(C_yearly[YRS_spinup:], 'C', paramind, cellind, version, shuffle, parameter)
#            fig_ts(N_yearly[YRS_spinup:], 'N', paramind, cellind, version, shuffle, parameter)
            #fig_ts(est_yearly[YRS_spinup:], 'est', paramind, cellind, version, shuffle, parameter)
        #    fig_ts(fpc_yearly[YRS_spinup:], 'fpc', paramind, cellind, version, shuffle, parameter)
            #fig_ts(mort_yearly[YRS_spinup:], 'mort', paramind, cellind, version, shuffle, parameter)


### N scatter
#            fig_N_scatter(N_beforeest_yearly[YRS_spinup:],est_yearly[YRS_spinup:],'N','E') #,paramind, cellind, version, shuffle, parameter)
#            fig_N_scatter(N_beforeest_yearly[YRS_spinup:],est_yearly[YRS_spinup:]+Nmort_yearly[YRS_spinup:]+adjN_yearly[YRS_spinup],'N','E+M+A') #,paramind, cellind, version, shuffle, parameter)
#        
#        #    fig_N_scatter(fpc_beforeest_yearly[YRS_spinup:],est_yearly[YRS_spinup:],'fpc','est',paramind, cellind, version, shuffle, parameter)
#        #    fig_N_scatter(N_beforeest_yearly[YRS_spinup:],fpc_beforeest_yearly[YRS_spinup:],'N','fpc',paramind, cellind, version, shuffle, parameter)
#      #      fig_N_scatter(N_beforeest_yearly[YRS_spinup:],Nmort_yearly[YRS_spinup:],'N','M',paramind, cellind, version, shuffle, parameter)
#      #      fig_N_scatter(N_beforeest_yearly[YRS_spinup:],adjN_yearly[YRS_spinup:],'N','A',paramind, cellind, version, shuffle, parameter)
#
#         #   fig_N_scatter(N_beforeest_yearly[YRS_spinup:],Nmort_yearly[YRS_spinup:],'N','Nmort',paramind, cellind, version, shuffle, parameter)
#        
#           
            ### Allocation scatter
            #vary_list={'AS',} # 'ASH', 'ALR', 'ASHfrac',} # 'ALRfrac',
            #varx_list={'A',}
            #varz='drought1'

            
            ######## allocation scatter
#            #for varx in varx_list:
##            #    for vary in vary_list:
#            fig_A_scatter(A_yearly[YRS_spinup:],AL_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','AL','drought')
#            fig_A_scatter(A_yearly[YRS_spinup:],AR_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','AR','drought')
            fig_A_scatter(A_yearly[YRS_spinup:],AS_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','AS','drought')
#            fig_A_scatter(A_yearly[YRS_spinup:],AH_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','AH','drought')

            fig_A_scatter(A_yearly[YRS_spinup:],AL_yearly[YRS_spinup:]+AR_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','ALR','drought')            
            fig_A_scatter(A_yearly[YRS_spinup:],AS_yearly[YRS_spinup:]+AH_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','ASH','drought')
            fig_A_scatter(A_yearly[YRS_spinup:],(AS_yearly[YRS_spinup:]+AH_yearly[YRS_spinup:])/A_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','ASHfrac','drought')
            
#            fig_A_scatter(Aperm2_yearly[YRS_spinup:],ALperm2_yearly[YRS_spinup:]+ARperm2_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','ALR','drought')
#            fig_A_scatter(Aperm2_yearly[YRS_spinup:],ASperm2_yearly[YRS_spinup:]+AHperm2_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','ASH','drought')
#            fig_A_scatter(Aperm2_yearly[YRS_spinup:],(ASperm2_yearly[YRS_spinup:]+AHperm2_yearly[YRS_spinup:])/Aperm2_yearly[YRS_spinup:],drought_yearly[YRS_spinup:], 'A','ASHfrac','drought')     

            
##          ### eval (compare to LPJ):
#            fig_ts_eval(AL_yearly[YRS_spinup:], AL_LPJ[YRS_spinup:, cellind], 'AL', paramind, cellind, version, shuffle, parameter)
#            fig_ts_eval(Cperm2_yearly[YRS_spinup:], Cperm2_LPJ[YRS_spinup:, cellind], 'Cperm2', paramind, cellind, version, shuffle, parameter)
#            fig_ts_eval(N_yearly[YRS_spinup:], N_afteryear_LPJ[YRS_spinup:, cellind], 'N', paramind, cellind, version, shuffle, parameter)
##            fig_ts_eval(est_yearly[YRS_spinup:], est_LPJ[YRS_spinup:, cellind], 'est', paramind, cellind, version, shuffle, parameter)
##            fig_ts_eval(mort_yearly[YRS_spinup:], mort_LPJ[YRS_spinup:, cellind], 'mort', paramind, cellind, version, shuffle, parameter)
##            fig_ts_eval(fpc_yearly[YRS_spinup:], fpc_LPJ[YRS_spinup:, cellind], 'fpc', paramind, cellind, version, shuffle, parameter)
        