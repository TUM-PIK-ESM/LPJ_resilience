#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 17:05:13 2023

@author: bathiany
"""
import numpy as np
import random
import netCDF4 as nc
import sys
sys.path.insert(1, '/home/bathiany/Projects/Vegetation_resilience_indicators/reduced_models/Cbalance/functions')
import LPJ_functions
import par
from ACF import ACF_fun
from ACF import ts_separate
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})
np.seterr(invalid='ignore')

LPJdatapath='/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/'


#models=("A2N1","A1N1",)
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

#shuffle=1

### paper
#models=("A1N1",)
#variants=("M1E1A1",)
#shuffle=1
##parameters=("Nvarfrac","varfrac")
#parameters=("Nvarfrac",)


## show that varfrac even reverses pattern when using fixed N
#models=("A1N3",)
#variants=("M0E1A1",)  # run again with 52000 years
#shuffle=1
#parameters=("varfrac",)

## for all versions plot
#versions=("A1N1M1E1A1S1","A1N1M0E1A1S1","A1N3M0E1A1S1","A2N1M0E1A1S1")
#models=("A1N1","A1N3","A2N1")
#variants=("M0E1A1",)
shuffle=1
parameters=("Nvarfrac",)

models=("A1N1",)
variants=("M0E1A1",)


### A0: prescribe allocated fluxes: models=("A0N1",)
#models=("A0N0","A0N1","A0N3")
#variants=("M1E1A1","M0E1A1",)
#shuffle=1
#parameters=("varfrac",)
#
#models=("A0N3",)
#variants=("M0E1A1",)
#shuffle=1
#parameters=("varfrac",)


##parameters=("ASHfrac")


########

dpi=500
exp="sba0084"      # exp from which data is read in
save_output=1

PFT=1
index_PS=0

Nparams=7
par.set_fixed_parameters(PFT)

PFTind=PFT-1

Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]
Pvec=np.array(Pvec)                 # in mm/year

#Nprec=19   # 
cellind_ini=0   # remove P=0, starting from 1 (200 mm/yr, i.e. cellind=1)
cellind_fin=20
Pvec=Pvec[cellind_ini:cellind_fin] 

#Pvec=Pvec[1:Nprec+1] # remove P=0, and cut off values above 2000mm; starting from 1 (200 mm/yr, i.e. cellind=1)

Nprec=cellind_fin-cellind_ini


YRS_spinup=2000      

if shuffle==0 or shuffle==3:
    YRS=10000
else:
    YRS=52000         #  works thanks to shuffling

#####
YRS_LPJ=YRS
if (YRS>10000):
    YRS_LPJ=10000

### test: short time periods only
#YRS=4000
            

########## list of output diagnostics
diagvar_list=("Cperm2","C","N","AC10_Cperm2","height",\
              "AC30_Cperm2","AC10_C","AC30_C","AC10_LSperm2",\
              "AC30_LSperm2","AC1trans_Cperm2","AC30trans_Cperm2",\
              "AC3trans_Cperm2","AC5trans_Cperm2","Tau",\
              "adjNrel","drought1","drought2","drought","fpcstd",\
              "mort","turnover","bm_delta","est","adjcount","estcount","mortest",\
              "estvsN","estvsC","fpc","Hfrac","AHfrac","SHfrac",\
              "ASHfrac","ACFmin","asym","AC10_N","Cperm2std",\
              "Cstd","Nstd", "eststd","NstdperN","AC10_est",\
              "estmin","estmax","VarN", "VarC",\
              "corr_Nfpc","corr_Nest","corr_fpcest",\
              "corr_CN_lag0", "cov_CN_lag0",\
              "corr_CN_lag10", "corr_NC_lag10",\
              "cov_CN_lag10", "cov_NC_lag10",\
              "BGZ1_lag10","BGZ2_lag10",\
              "BGZ3_lag10","BGZ4_lag10",\
              "BGZ5_lag10","BGZ6_lag10",\
              "BGZ7_lag10","BGZ8_lag10",\
              "BGZ9_lag10","BGZ10","BGZ_lag10",\
              "BGN1","BGN2","BGN3","BGN4","BGN5","BGN6","BGN7","BGN",\
              "AC10_Cperm2_BG","cov_Cperm2_lag10","VarCperm2","lambdaN")


## correlations and covariances in case they are needed for diagnostics:
## VarN, VarC, N2, C2, 
# corr_CN_lag0, cov_CN_lag0,
# corr_CN_lag10, corr_NC_lag10, corr_CN_lag30, corr_NC_lag30,
# cov_CN_lag10,  cov_NC_lag10,  cov_CN_lag30,  cov_NC_lag30,


## Bohrnstedt&Goldberger terms, numerator (Z) 10 Terms (lag-dependent except no 10), BGZi
# and denumerator (N), 7 Terms, BGNi:
#"BGZ1_lag10","BGZ1_lag30","BGZ2_lag10","BGZ2_lag30",
#"BGZ3_lag10","BGZ3_lag30","BGZ4_lag10","BGZ4_lag30",
#"BGZ5_lag10","BGZ5_lag30","BGZ6_lag10","BGZ6_lag30",
#"BGZ7_lag10","BGZ7_lag30","BGZ8_lag10","BGZ8_lag30",
#"BGZ9_lag10","BGZ9_lag30","BGZ10_lag10","BGZ10_lag30",
#"BGN1","BGN2","BGN3","BGN4","BGN5","BGN6","BGN7")

# for evaluation, the total result (which is identical to AC10 as it should be):
#AC10_Cperm2_BG_diag

############ unused correlations:
#               "corr_CN_P_diag", "corr_Cperm2N_P_diag", "corr_Cperm2C_P_diag",\
#               "corr_CN_S_diag", "corr_Cperm2N_S_diag", "corr_Cperm2C_S_diag",\
#               "corr_CN_K_diag", "corr_Cperm2N_K_diag", "corr_Cperm2C_K_diag")

##########


## loop over model versions
for model in models:
    for variant in variants:
        version=model+variant+'S'+str(shuffle)
        adj_switch=int(variant[5])      
        
        ########## read driver data from full LPJ model
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
    
    
        ## allocation fluxes
        # this is per individual, will be translated to per m2 using the fixed prescribed N
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
        #N_afteryear_LPJ_mean=np.mean(N_afteryear_LPJ, axis=0)
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
        NAC_LPJ=np.zeros(np.shape(N_afteryear_LPJ)[1])
        # NAC at lag 10 shall be = NAC1**10
        for ind in range (0,np.shape(N_afteryear_LPJ)[1]):
        #    NAC10_LPJ=ACF_fun(N_afteryear_LPJ[:,ind], 11)[10]
        #    NAC_LPJ[ind]=np.power(NAC10_LPJ,1/10)
            NAC30_LPJ=ACF_fun(N_afteryear_LPJ[:,ind], 31)[30]
            NAC_LPJ[ind]=np.power(NAC30_LPJ,1/30)
        sigma_noise=Nstd_LPJ*np.sqrt(1-NAC_LPJ)
    

        for parameter in parameters:
            #############################
    
    
    
            ####### initialise matrices for output diagnostics
            AC10_Cperm2_diag=np.zeros((Nparams,Nprec))
            AC30_Cperm2_diag=np.zeros((Nparams,Nprec))
            AC10_C_diag=np.zeros((Nparams,Nprec))
            AC30_C_diag=np.zeros((Nparams,Nprec))
            
            AC10_LSperm2_diag=np.zeros((Nparams,Nprec))
            AC30_LSperm2_diag=np.zeros((Nparams,Nprec))
            
            AC1trans_Cperm2_diag=np.zeros((Nparams,Nprec))
            AC2trans_Cperm2_diag=np.zeros((Nparams,Nprec))
            AC3trans_Cperm2_diag=np.zeros((Nparams,Nprec))
            AC5trans_Cperm2_diag=np.zeros((Nparams,Nprec))
            AC10trans_Cperm2_diag=np.zeros((Nparams,Nprec))
            AC30trans_Cperm2_diag=np.zeros((Nparams,Nprec))
    
    
            Tau_diag=np.zeros((Nparams,Nprec))
    
            L_diag=np.zeros((Nparams,Nprec))
            R_diag=np.zeros((Nparams,Nprec))
            S_diag=np.zeros((Nparams,Nprec))
            H_diag=np.zeros((Nparams,Nprec))
            Cperm2_diag=np.zeros((Nparams,Nprec))
            C_diag=np.zeros((Nparams,Nprec))
            Czero_diag=np.zeros((Nparams,Nprec))
            height_diag=np.zeros((Nparams,Nprec))
            
            Hfrac_diag=np.zeros((Nparams,Nprec))
            AHfrac_diag=np.zeros((Nparams,Nprec))
            SHfrac_diag=np.zeros((Nparams,Nprec))
            ASHfrac_diag=np.zeros((Nparams,Nprec))
            mort_diag=np.zeros((Nparams,Nprec))
            bm_delta_diag=np.zeros((Nparams,Nprec))
            turnover_diag=np.zeros((Nparams,Nprec))
            est_diag=np.zeros((Nparams,Nprec))
            adjN_diag=np.zeros((Nparams,Nprec))
            adjNrel_diag=np.zeros((Nparams,Nprec))
            drought1_diag=np.zeros((Nparams,Nprec))
            drought2_diag=np.zeros((Nparams,Nprec))
    
            adjcount_diag=np.zeros((Nparams,Nprec))
            estcount_diag=np.zeros((Nparams,Nprec))
            
            mortest_diag=np.zeros((Nparams,Nprec))
            
            estvsN_diag=np.zeros((Nparams,Nprec))
            estvsC_diag=np.zeros((Nparams,Nprec))
    
            fpc_diag=np.zeros((Nparams,Nprec))
            fpcstd_diag=np.zeros((Nparams,Nprec))
            
            N_diag=np.zeros((Nparams,Nprec))
            Nstd_diag=np.zeros((Nparams,Nprec))
            eststd_diag=np.zeros((Nparams,Nprec))
            estmin_diag=np.zeros((Nparams,Nprec))
            estmax_diag=np.zeros((Nparams,Nprec))
    
            NstdperN_diag=np.zeros((Nparams,Nprec))
            AC10_N_diag=np.zeros((Nparams,Nprec))
    
            AC10_est_diag=np.zeros((Nparams,Nprec))
            AC30_est_diag=np.zeros((Nparams,Nprec))
    
            lambdaN_diag=np.zeros((Nparams,Nprec))
    
            asym_diag=np.zeros((Nparams,Nprec))
            ACFmin_diag=np.zeros((Nparams,Nprec))
            
            Cstd_diag=np.zeros((Nparams,Nprec))
            Cperm2std_diag=np.zeros((Nparams,Nprec))
    
            corr_Nfpc_diag=np.zeros((Nparams,Nprec))
            corr_Nest_diag=np.zeros((Nparams,Nprec))
            corr_fpcest_diag=np.zeros((Nparams,Nprec))
            
    #        corr_CN_P_diag=np.zeros((Nparams,Nprec))
    #        corr_Cperm2N_P_diag=np.zeros((Nparams,Nprec))
    #        corr_Cperm2C_P_diag=np.zeros((Nparams,Nprec))
    #        corr_CN_S_diag=np.zeros((Nparams,Nprec))
    #        corr_Cperm2N_S_diag=np.zeros((Nparams,Nprec))
    #        corr_Cperm2C_S_diag=np.zeros((Nparams,Nprec))
    #        corr_CN_K_diag=np.zeros((Nparams,Nprec))
    #        corr_Cperm2N_K_diag=np.zeros((Nparams,Nprec))
    #        corr_Cperm2C_K_diag=np.zeros((Nparams,Nprec))
    #        
    
            VarN_diag=np.zeros((Nparams,Nprec))
            VarC_diag=np.zeros((Nparams,Nprec))  
            #N2_diag=np.zeros((Nparams,Nprec))  
            #C2_diag=np.zeros((Nparams,Nprec)) 
            corr_CN_lag0_diag=np.zeros((Nparams,Nprec))  
            cov_CN_lag0_diag=np.zeros((Nparams,Nprec)) 
            corr_CN_lag10_diag=np.zeros((Nparams,Nprec))  
            corr_NC_lag10_diag=np.zeros((Nparams,Nprec))  
            #corr_CN_lag30_diag=np.zeros((Nparams,Nprec))  
            #corr_NC_lag30_diag=np.zeros((Nparams,Nprec)) 
            cov_CN_lag10_diag=np.zeros((Nparams,Nprec))  
            cov_NC_lag10_diag=np.zeros((Nparams,Nprec))  
            #cov_CN_lag30_diag=np.zeros((Nparams,Nprec))  
            #cov_NC_lag30_diag=np.zeros((Nparams,Nprec)) 
            BGZ1_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ1_lag30_diag=np.zeros((Nparams,Nprec))  
            BGZ2_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ2_lag30_diag=np.zeros((Nparams,Nprec)) 
            BGZ3_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ3_lag30_diag=np.zeros((Nparams,Nprec))  
            BGZ4_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ4_lag30_diag=np.zeros((Nparams,Nprec)) 
            BGZ5_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ5_lag30_diag=np.zeros((Nparams,Nprec))  
            BGZ6_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ6_lag30_diag=np.zeros((Nparams,Nprec)) 
            BGZ7_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ7_lag30_diag=np.zeros((Nparams,Nprec))  
            BGZ8_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ8_lag30_diag=np.zeros((Nparams,Nprec)) 
            BGZ9_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ9_lag30_diag=np.zeros((Nparams,Nprec))  
            BGZ10_diag=np.zeros((Nparams,Nprec))  
            BGZ_lag10_diag=np.zeros((Nparams,Nprec))  
            #BGZ_lag30_diag=np.zeros((Nparams,Nprec))  
    
            BGN1_diag=np.zeros((Nparams,Nprec))  
            BGN2_diag=np.zeros((Nparams,Nprec))  
            BGN3_diag=np.zeros((Nparams,Nprec))  
            BGN4_diag=np.zeros((Nparams,Nprec))  
            BGN5_diag=np.zeros((Nparams,Nprec))  
            BGN6_diag=np.zeros((Nparams,Nprec))  
            BGN7_diag=np.zeros((Nparams,Nprec)) 
            BGN_diag=np.zeros((Nparams,Nprec))  
    
            AC10_Cperm2_BG_diag=np.zeros((Nparams,Nprec)) 
    
           
            cov_Cperm2_lag10_diag=np.zeros((Nparams,Nprec))  
            VarCperm2_diag=np.zeros((Nparams,Nprec))  
           
                    
            
            ## initial condition from LPJ, wet cell:
            L, R, S, H, D, N, fpc, height = 9057.575, 9157.507, 107325.516, 140601.61, 0.0, 0.06495886, 0.95, 17.915508
            
            for paramind in range(0,Nparams):
            
                ## set parameter values
                par.set_variable_parameters(Nparams, parameter, paramind)
                Nvarfrac=par.Nvarfrac
                varfrac=par.varfrac
                estfrac=par.estfrac
                AHfrac=par.AHfrac
                ASHfrac=par.ASHfrac
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
                    # ALfrac and ARfrac are unchanged, but more AH reduces AS.                                        
                    ALfrac=ALfrac_LPJ
                    ARfrac=ARfrac_LPJ
                    ASfrac=1-ALfrac-ARfrac-AHfrac
    
                elif parameter=="ASHfrac": 
                    # ASHfrac is prescribed, the rest is scaled based on average fracs.
                    # ASHfrac is what's shown in Allocation scatter plot, panel d, horizonal dashed red line.
                    
                    LofLR=ALfrac_LPJ/(ALfrac_LPJ+ARfrac_LPJ)
                    SofSH=ASfrac_LPJ/(ASfrac_LPJ+AHfrac_LPJ)
                    ALRfrac=1-ASHfrac
                    ALfrac=ALRfrac*LofLR
                    ARfrac=ALRfrac*(1-LofLR)
                    ASfrac=ASHfrac*SofSH
                    AHfrac=ASHfrac*(1-SofSH)
 
                    #print(LofLR, SofSH, ALRfrac, ALfrac, ARfrac, ASfrac, AHfrac)
                   
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
                    N_afteryear_input=Nvarfrac*(N_afteryear_LPJ-Nmean_LPJ)+Nmean_LPJ
                    N_beforeest_input=Nvarfrac*(N_beforeest_LPJ-N_beforeest_LPJ_mean)+N_beforeest_LPJ_mean
                    N_afterest_input=Nvarfrac*(N_afterest_LPJ-N_afterest_LPJ_mean)+N_afterest_LPJ_mean
                    N_beforemort_input=Nvarfrac*(N_beforemort_LPJ-N_beforemort_LPJ_mean)+N_beforemort_LPJ_mean
                
                
                
                for cellind in range(cellind_ini,cellind_fin):
    
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
    
                    #print(est_prescribed)
                    
                    # diagnosis           
                    mort_yearly=np.zeros(YRS)
                    bm_delta_yearly=np.zeros(YRS)
                    turnover_yearly=np.zeros(YRS)
                    est_yearly=np.zeros(YRS)
                    adjcount_yearly=np.zeros(YRS)
                    estcount_yearly=np.zeros(YRS)
                    fpc_beforeest_yearly=np.zeros(YRS)
                    N_beforeest_yearly=np.zeros(YRS)
        
                    adjN_yearly=np.zeros(YRS)
                    adjNrel_yearly=np.zeros(YRS)
                    drought1_yearly=np.zeros(YRS)
                    drought2_yearly=np.zeros(YRS)
                    C_yearly=np.zeros(YRS)
                    Cperm2_yearly=np.zeros(YRS)
                    LSperm2_yearly=np.zeros(YRS)
    
                    height_yearly=np.zeros(YRS)
        
                    estvsN_yearly=np.zeros(YRS)
                    estvsC_yearly=np.zeros(YRS)
    
    
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
            
                    Czero_yearly=np.zeros(YRS)
            
                    # inicond in case last cell ended up in a mess
                    L, R, S, H, D, N, fpc, height = 9057.575, 9157.507, 107325.516, 140601.61, 0.0, 0.06495886, 0.95, 17.915508
                    
                    for year in range(1,YRS):
            
                        
                        
                    ####### take care of shuffling type (randomisation of inputs)
                        
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
                            
                            
                        ## set diagnostic outputs to 0, then models will overwrite those that they compute
                        drought1=0
                        drought2=0
                        adjN=0
                        adjNrel=0
                        mort=0
                        bm_delta=0
                        turnover=0                    
                        est=0
                        estvsN=0
                        estvsC=0
                        adj_count=0
                        est_count=0
                        fpc_beforeest=0
                        N_beforeest=0                    
    
    
                        ### compute outputs that are not directly generated by the model in some variants
                    
                        if (version[0:2]=="A0"):        # or version[0:2]=="A2"                   
                            ALperm2=AL_input[year_A, cellind]*N_afteryear_LPJ[year_N-1, cellind]
                            ARperm2=AR_input[year_A, cellind]*N_afteryear_LPJ[year_N-1, cellind]
                            ASperm2=AS_input[year_A, cellind]*N_afteryear_LPJ[year_N-1, cellind]
                            AHperm2=AH_input[year_A, cellind]*N_afteryear_LPJ[year_N-1, cellind]                 
            
                        if (version[0:2]=="A2" or version[0:2]=="A3"):
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
                        elif (version[2:4]=="N3"):   # constant N over time
                            N=Nmean_LPJ[cellind]
                        elif (version[2:4]=="N7"):   # AR1 process
                            if shuffle==0:
                                N = NAC_LPJ[cellind]*(N-Nmean_LPJ[cellind]) + random.gauss(Nmean_LPJ[cellind],sigma_noise[cellind])
                            else:
                                N = random.gauss(Nmean_LPJ[cellind],Nstd_LPJ[cellind])
                
                        ### take care of prescribed mortaility and/or establishment
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
                        elif (version[6:8]=="E0"):   # fixed at global constant
                            est_prescribed=0.01
                        elif version[6:8]=="E1":     # interactive
                            est_prescribed=-999
                        elif (version[6:8]=="E2"):   # fixed to local mean
                            est_prescribed=np.mean(est_LPJ[:,cellind])
                        elif (version[6:8]=="E3"):   # shuffled in time (destroys correlations)
                            est_prescribed=est_LPJ[year_E,cellind]                

                        ## two options, which are theoretically identical if input is not shuffled, but which matter for numerical stability:
                        bm_forcing = bm_inc_ind_input[year_A, cellind]  # this is the default mode; passes bm_inc per individual from LPJ
                        #bm_forcing = bm_inc_input[year_A, cellind]/N  # alternative mode; passes bm_inc per m2 from LPJ

                        ############################## run the model
                        if model=='A1N1':
                            L, R, S, H, D, N, fpc, height, AL, AR, AS, AH, drought1, drought2, \
                                adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest = \
                                   LPJ_functions.LPJ_Cbalance_A1N1(L, R, S, H, \
                                      D, N, height, bm_forcing, \
                                      wscal_input[year_A, cellind], fL, fR, fS, Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch)
                            ## in case E and M are fixed, there is no variability in N. In this case, add noise:
                            if (parameter=="est" or version[6:8]=="E0" or version[6:8]=="E2") and (parameter=="mort" or version[4:6]=="M0" or version[4:6]=="M2"):
                                N=N+random.gauss(0,Nstd_LPJ[cellind])  
                                
                        elif model=='A0N1':  
                            L, R, S, H, N, fpc, height, \
                                adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest = \
                                LPJ_functions.LPJ_Cbalance_A0N1(L, R, S, H, \
                                   N, height, AL, AR, AS, AH, fL, fR, fS, Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch)
        
                            ## only for diagnostic purpose:                            
                            AL=ALperm2/N
                            AR=ARperm2/N
                            AS=ASperm2/N
                            AH=AHperm2/N
                                            
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
        
    
                        ############## save this year's output
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
            
                        N_yearly[year]=N
                        fpc_yearly[year]=fpc
    
                        height_yearly[year]=height
                        
                        drought1_yearly[year]=drought1
                        drought2_yearly[year]=drought2
                        adjcount_yearly[year]=adj_count
                        estcount_yearly[year]=est_count
    
                        fpc_beforeest_yearly[year]=fpc_beforeest
                        N_beforeest_yearly[year]=N_beforeest
                        
                        adjN_yearly[year]=adjN
                        adjNrel_yearly[year]=adjNrel
                        mort_yearly[year]=mort
                        bm_delta_yearly[year]=bm_delta
                        turnover_yearly[year]=turnover
                        est_yearly[year]=est
                        estvsN_yearly[year]=estvsN
    
                        estvsC_yearly[year]=par.Csapl/(L+R+S+H)*estvsN*estfrac
    
            
                        if L<1e-9 or R<1e-9 or S<1e-9 or H<1e-9:
                            Czero_yearly[year]=1
            
                    #### end of time loop
                    
                    
                    Cperm2_anomaly_bw40, Cperm2_smoothed_bw40 = ts_separate(Cperm2_yearly[YRS_spinup:], bw=40)
                    Cperm2_anomaly_bw100, Cperm2_smoothed_bw100 = ts_separate(Cperm2_yearly[YRS_spinup:], bw=100)
                            
                    
                    
                    ############################## write out the diagnostic results
                    ## ACFs
                    ACF_Cperm2=ACF_fun(Cperm2_yearly[YRS_spinup:], 31)
                    AC10_Cperm2_diag[paramind,cellind]=ACF_Cperm2[10]
                    AC30_Cperm2_diag[paramind,cellind]=ACF_Cperm2[30]
            
                    ACF_C=ACF_fun(C_yearly[YRS_spinup:], 31)
                    AC10_C_diag[paramind,cellind]=ACF_C[10]
                    AC30_C_diag[paramind,cellind]=ACF_C[30]
            
                    ## just leaf and sapwood because satellites (VOD) would not see
                    # much H and R
                    ACF_LSperm2=ACF_fun((LSperm2_yearly[YRS_spinup:]), 31)
                    AC10_LSperm2_diag[paramind,cellind]=ACF_LSperm2[10]
                    AC30_LSperm2_diag[paramind,cellind]=ACF_LSperm2[30]
            
                    # high-freq part:
                    ACF_Cperm2_bw40=ACF_fun(Cperm2_anomaly_bw40, 11)
                    ACF_Cperm2_bw100=ACF_fun(Cperm2_anomaly_bw100, 31)
    
                    AC1trans_Cperm2_diag[paramind,cellind]=ACF_Cperm2_bw40[1]
                    AC3trans_Cperm2_diag[paramind,cellind]=ACF_Cperm2_bw40[3]
                    AC5trans_Cperm2_diag[paramind,cellind]=ACF_Cperm2_bw40[5]
                    AC10trans_Cperm2_diag[paramind,cellind]=ACF_Cperm2_bw40[10]
                    AC30trans_Cperm2_diag[paramind,cellind]=ACF_Cperm2_bw100[30]
                
                    #Tau_diag[paramind,cellind]=np.mean(C_yearly[YRS_spinup:])/np.mean(bm_inc_ind_input[YRS_spinup:, cellind])
        
                    L_diag[paramind,cellind]=np.mean(L_yearly[YRS_spinup:])
                    R_diag[paramind,cellind]=np.mean(R_yearly[YRS_spinup:])
                    S_diag[paramind,cellind]=np.mean(S_yearly[YRS_spinup:])
                    H_diag[paramind,cellind]=np.mean(H_yearly[YRS_spinup:])
                    Cperm2_diag[paramind,cellind]=np.mean(Cperm2_yearly[YRS_spinup:])
                    C_diag[paramind,cellind]=np.mean(C_yearly[YRS_spinup:])
                    height_diag[paramind,cellind]=np.mean(height_yearly[YRS_spinup:])
    
                    Czero_diag[paramind,cellind]=np.mean(Czero_yearly[YRS_spinup:])
                
                    Cstd_diag[paramind,cellind]=np.std(C_yearly[YRS_spinup:])
                    Cperm2std_diag[paramind,cellind]=np.std(Cperm2_yearly[YRS_spinup:])
                    N_diag[paramind,cellind]=np.mean(N_yearly[YRS_spinup:])
                    Nstd_diag[paramind,cellind]=np.std(N_yearly[YRS_spinup:])
                    eststd_diag[paramind,cellind]=np.std(est_yearly[YRS_spinup:])
    
                    estmin_diag[paramind,cellind]=np.min(est_yearly[YRS_spinup:])
                    estmax_diag[paramind,cellind]=np.max(est_yearly[YRS_spinup:])
    
                    corr_Nfpc_diag[paramind,cellind]=np.corrcoef(N_beforeest_yearly[YRS_spinup:], fpc_beforeest_yearly[YRS_spinup:])[0,1]
                    corr_Nest_diag[paramind,cellind]=np.corrcoef(N_beforeest_yearly[YRS_spinup:], est_yearly[YRS_spinup:])[0,1]
                    corr_fpcest_diag[paramind,cellind]=np.corrcoef(fpc_beforeest_yearly[YRS_spinup:], est_yearly[YRS_spinup:])[0,1]
    
                    NstdperN_diag[paramind,cellind]=np.std(N_yearly[YRS_spinup:])/np.mean(N_yearly[YRS_spinup:])
    
                    ACF_N=ACF_fun(N_yearly[YRS_spinup:], 1000)
                    AC10_N_diag[paramind,cellind]=ACF_N[10]
    
                    delta_lower=np.median(N_yearly[YRS_spinup:])-np.min(N_yearly[YRS_spinup:])
                    delta_upper=np.max(N_yearly[YRS_spinup:])-np.median(N_yearly[YRS_spinup:])
                    if delta_lower>delta_upper:                        
                        asym_diag[paramind,cellind]=delta_lower/delta_upper
                    else:
                        asym_diag[paramind,cellind]=delta_upper/delta_lower
                    
                    ACFmin_diag[paramind,cellind]=np.min(ACF_N[:])
            
                    SHfrac_yearly=(S_yearly+H_yearly)/(L_yearly+R_yearly+S_yearly+H_yearly)
                    ASHfrac_yearly=(AS_yearly+AH_yearly)/(AL_yearly+AR_yearly+AS_yearly+AH_yearly)
            
                    Hfrac_yearly=H_yearly/(L_yearly+R_yearly+S_yearly+H_yearly)
                    AHfrac_yearly=AH_yearly/(AL_yearly+AR_yearly+AS_yearly+AH_yearly)
            
                    Hfrac_diag[paramind,cellind]=np.nanmean(Hfrac_yearly[YRS_spinup:])
                    AHfrac_diag[paramind,cellind]=np.nanmean(AHfrac_yearly[YRS_spinup:])
            
                    SHfrac_diag[paramind,cellind]=np.nanmean(SHfrac_yearly[YRS_spinup:])
                    ASHfrac_diag[paramind,cellind]=np.nanmean(ASHfrac_yearly[YRS_spinup:])
                    
                    mort_diag[paramind,cellind]=np.mean(mort_yearly[YRS_spinup:])
                    bm_delta_diag[paramind,cellind]=np.mean(bm_delta_yearly[YRS_spinup:])
                    turnover_diag[paramind,cellind]=np.mean(turnover_yearly[YRS_spinup:])
                    est_diag[paramind,cellind]=np.mean(est_yearly[YRS_spinup:])
                    mortest_diag[paramind,cellind]=mort_diag[paramind,cellind]+est_diag[paramind,cellind]/N_diag[paramind,cellind]
    
                    adjN_diag[paramind,cellind]=-np.mean(adjN_yearly[YRS_spinup:])
                    adjNrel_diag[paramind,cellind]=-np.mean(adjNrel_yearly[YRS_spinup:])
                    drought1_diag[paramind,cellind]=np.mean(drought1_yearly[YRS_spinup:])
                    drought2_diag[paramind,cellind]=np.mean(drought2_yearly[YRS_spinup:])
    
                    estvsN_diag[paramind,cellind]=np.mean(estvsN_yearly[YRS_spinup:])
                    estvsC_diag[paramind,cellind]=np.mean(estvsC_yearly[YRS_spinup:])
    
                    estcount_diag[paramind,cellind]=np.mean(estcount_yearly[YRS_spinup:])
                    adjcount_diag[paramind,cellind]=np.mean(adjcount_yearly[YRS_spinup:])
    
    
                    fpc_diag[paramind,cellind]=np.mean(fpc_yearly[YRS_spinup:])
                    fpcstd_diag[paramind,cellind]=np.std(fpc_yearly[YRS_spinup:])
    
                    ACF_est=ACF_fun(est_yearly[YRS_spinup:], 31)
                    AC10_est_diag[paramind,cellind]=ACF_est[10]
                    #AC30_est_diag[paramind,cellind]=ACF_est[30]
    
    
                    DNmean=np.mean(N_yearly[YRS_spinup+1:]-N_yearly[YRS_spinup:-1])
                    lambdaN_diag[paramind,cellind]=DNmean/N_diag[paramind,cellind]  # should be ~ mort+est/N
    
    
                    #print(paramind,cellind,AC30_diag[paramind,cellind])
    #
    #                corr_CN_P_diag[paramind,cellind]=scipy.stats.pearsonr(C_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0]
    #                corr_Cperm2N_P_diag[paramind,cellind]=scipy.stats.pearsonr(Cperm2_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0]
    #                corr_Cperm2C_P_diag[paramind,cellind]=scipy.stats.pearsonr(Cperm2_yearly[YRS_spinup:], C_yearly[YRS_spinup:])[0]
    
    #                corr_CN_S_diag[paramind,cellind]=scipy.stats.spearmanr(C_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0]
    #                corr_Cperm2N_S_diag[paramind,cellind]=scipy.stats.spearmanr(Cperm2_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0]
    #                corr_Cperm2C_S_diag[paramind,cellind]=scipy.stats.spearmanr(Cperm2_yearly[YRS_spinup:], C_yearly[YRS_spinup:])[0]
    
    #                corr_CN_K_diag[paramind,cellind]=scipy.stats.kendalltau(C_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0]
    #                corr_Cperm2N_K_diag[paramind,cellind]=scipy.stats.kendalltau(Cperm2_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0]
    #                corr_Cperm2C_K_diag[paramind,cellind]=scipy.stats.kendalltau(Cperm2_yearly[YRS_spinup:], C_yearly[YRS_spinup:])[0]
    
    
    
    
                    ###### diagnostics related to indiv terms determining the AC of Cperm2
    
                    Nanom=N_yearly-np.mean(N_yearly[YRS_spinup:])
                    Canom=C_yearly-np.mean(C_yearly[YRS_spinup:])
            
                    VarN_diag[paramind,cellind]=np.var(N_yearly[YRS_spinup:])
                    VarC_diag[paramind,cellind]=np.var(C_yearly[YRS_spinup:])
                    VarCperm2_diag[paramind,cellind]=np.var(Cperm2_yearly[YRS_spinup:])
    
                    #N2_diag[paramind,cellind]=np.mean(N_yearly[YRS_spinup:]*N_yearly[YRS_spinup:])
                    #C2_diag[paramind,cellind]=np.mean(C_yearly[YRS_spinup:]*C_yearly[YRS_spinup:])
    
                    corr_CN_lag0_diag[paramind,cellind]=np.corrcoef(C_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0,1]
                    cov_CN_lag0_diag[paramind,cellind]=np.cov(C_yearly[YRS_spinup:], N_yearly[YRS_spinup:])[0,1]
    
                    corr_CN_lag10_diag[paramind,cellind]=np.corrcoef(C_yearly[YRS_spinup+10:], N_yearly[YRS_spinup:-10])[0,1]
                    corr_NC_lag10_diag[paramind,cellind]=np.corrcoef(N_yearly[YRS_spinup+10:], C_yearly[YRS_spinup:-10])[0,1]
                    cov_CN_lag10_diag[paramind,cellind]=np.cov(C_yearly[YRS_spinup+10:], N_yearly[YRS_spinup:-10])[0,1]
                    cov_NC_lag10_diag[paramind,cellind]=np.cov(N_yearly[YRS_spinup+10:], C_yearly[YRS_spinup:-10])[0,1]
    
    #                corr_CN_lag30_diag[paramind,cellind]=np.corrcoef(C_yearly[YRS_spinup+30:], N_yearly[YRS_spinup:-30])[0,1]
    #                corr_NC_lag30_diag[paramind,cellind]=np.corrcoef(N_yearly[YRS_spinup+30:], C_yearly[YRS_spinup:-30])[0,1]
    #                cov_CN_lag30_diag[paramind,cellind]=np.cov(C_yearly[YRS_spinup+30:], N_yearly[YRS_spinup:-30])[0,1]
    #                cov_NC_lag30_diag[paramind,cellind]=np.cov(N_yearly[YRS_spinup+30:], C_yearly[YRS_spinup:-30])[0,1]
    
                    cov_CC_lag10=np.cov(C_yearly[YRS_spinup:-10], C_yearly[YRS_spinup+10:])[0,1]
                    cov_NN_lag10=np.cov(N_yearly[YRS_spinup:-10], N_yearly[YRS_spinup+10:])[0,1]
                    #cov_CC_lag30=np.cov(C_yearly[YRS_spinup:-30], C_yearly[YRS_spinup+30:])[0,1]
                    #cov_NN_lag30=np.cov(N_yearly[YRS_spinup:-30], N_yearly[YRS_spinup+30:])[0,1]
    
                    cov_Cperm2_lag10_diag[paramind,cellind]=np.cov(Cperm2_yearly[YRS_spinup:-10], Cperm2_yearly[YRS_spinup+10:])[0,1]
    
                    BGZ1_lag10=C_diag[paramind,cellind]*C_diag[paramind,cellind]*cov_NN_lag10
                    #BGZ1_lag30=C_diag[paramind,cellind]*C_diag[paramind,cellind]*cov_NN_lag30
                    BGZ2_lag10=C_diag[paramind,cellind]*N_diag[paramind,cellind]*cov_NC_lag10_diag[paramind,cellind]
                    #BGZ2_lag30=C_diag[paramind,cellind]*N_diag[paramind,cellind]*cov_NC_lag30_diag[paramind,cellind]
                    BGZ3_lag10=C_diag[paramind,cellind]*N_diag[paramind,cellind]*cov_CN_lag10_diag[paramind,cellind]
                    #BGZ3_lag30=C_diag[paramind,cellind]*N_diag[paramind,cellind]*cov_CN_lag30_diag[paramind,cellind]
                    BGZ4_lag10=N_diag[paramind,cellind]*N_diag[paramind,cellind]*cov_CC_lag10
                    #BGZ4_lag30=N_diag[paramind,cellind]*N_diag[paramind,cellind]*cov_CC_lag30
                    BGZ5_lag10=np.mean(Canom[YRS_spinup+10:]*Nanom[YRS_spinup+10:]*Canom[YRS_spinup:-10]*Nanom[YRS_spinup:-10])
                    #BGZ5_lag30=np.mean(Canom[YRS_spinup+30:]*Nanom[YRS_spinup+30:]*Canom[YRS_spinup:-30]*Nanom[YRS_spinup:-30])
                    BGZ6_lag10=C_diag[paramind,cellind]*np.mean(Nanom[YRS_spinup+10:]*Canom[YRS_spinup:-10]*Nanom[YRS_spinup:-10])
                    #BGZ6_lag30=C_diag[paramind,cellind]*np.mean(Nanom[YRS_spinup+30:]*Canom[YRS_spinup:-30]*Nanom[YRS_spinup:-30])
                    BGZ7_lag10=N_diag[paramind,cellind]*np.mean(Canom[YRS_spinup+10:]*Canom[YRS_spinup:-10]*Nanom[YRS_spinup:-10])
                    #BGZ7_lag30=N_diag[paramind,cellind]*np.mean(Canom[YRS_spinup+30:]*Canom[YRS_spinup:-30]*Nanom[YRS_spinup:-30])
                    BGZ8_lag10=C_diag[paramind,cellind]*np.mean(Canom[YRS_spinup+10:]*Nanom[YRS_spinup+10:]*Nanom[YRS_spinup:-10])
                    #BGZ8_lag30=C_diag[paramind,cellind]*np.mean(Canom[YRS_spinup+30:]*Nanom[YRS_spinup+30:]*Nanom[YRS_spinup:-30])
                    BGZ9_lag10=N_diag[paramind,cellind]*np.mean(Canom[YRS_spinup+10:]*Nanom[YRS_spinup+10:]*Canom[YRS_spinup:-10])
                    #BGZ9_lag30=N_diag[paramind,cellind]*np.mean(Canom[YRS_spinup+30:]*Nanom[YRS_spinup+30:]*Canom[YRS_spinup:-30])
                    BGZ10=(cov_CN_lag0_diag[paramind,cellind])**2

                    BGN1=C_diag[paramind,cellind]*C_diag[paramind,cellind]*VarN_diag[paramind,cellind]
                    BGN2=N_diag[paramind,cellind]*N_diag[paramind,cellind]*VarC_diag[paramind,cellind]
                    BGN3=np.mean(Nanom[YRS_spinup:]*Nanom[YRS_spinup:]*Canom[YRS_spinup:]*Canom[YRS_spinup:])
                    BGN4=2*C_diag[paramind,cellind]*np.mean(Canom[YRS_spinup:]*Nanom[YRS_spinup:]*Nanom[YRS_spinup:])
                    BGN5=2*N_diag[paramind,cellind]*np.mean(Canom[YRS_spinup:]*Canom[YRS_spinup:]*Nanom[YRS_spinup:])
                    BGN6=2*C_diag[paramind,cellind]*N_diag[paramind,cellind]*cov_CN_lag0_diag[paramind,cellind]
                    BGN7=-(cov_CN_lag0_diag[paramind,cellind])**2
    
                    BGZ_lag10_diag[paramind,cellind]=(BGZ1_lag10+BGZ2_lag10\
                                       +BGZ3_lag10+BGZ4_lag10\
                                       +BGZ5_lag10+BGZ6_lag10\
                                       +BGZ7_lag10+BGZ8_lag10\
                                       +BGZ9_lag10+BGZ10)
    
                    BGN_diag[paramind,cellind]=BGN1+BGN2\
                                         +BGN3+BGN4\
                                         +BGN5+BGN6+BGN7
    
                    AC10_Cperm2_BG_diag[paramind,cellind]=(BGZ1_lag10+BGZ2_lag10\
                                       +BGZ3_lag10+BGZ4_lag10\
                                       +BGZ5_lag10+BGZ6_lag10\
                                       +BGZ7_lag10+BGZ8_lag10\
                                       +BGZ9_lag10+BGZ10)\
                                       /(BGN1+BGN2\
                                         +BGN3+BGN4\
                                         +BGN5+BGN6+BGN7)
    #                AC30_Cperm2_BG_diag[paramind,cellind]=(BGZ1_lag30+BGZ2_lag30\
    #                                   +BGZ3_lag30+BGZ4_lag30\
    #                                   +BGZ5_lag30+BGZ6_lag30\
    #                                   +BGZ7_lag30+BGZ8_lag30\
    #                                   +BGZ9_lag30+BGZ10)\
    #                                   /(BGN1+BGN2\
    #                                     +BGN3+BGN4\
    #                                     +BGN5+BGN6+BGN7)
    
    
                    ### output of the terms, but normalised:
                    BGZ1_lag10_diag[paramind,cellind]=BGZ1_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ2_lag10_diag[paramind,cellind]=BGZ2_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ3_lag10_diag[paramind,cellind]=BGZ3_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ4_lag10_diag[paramind,cellind]=BGZ4_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ5_lag10_diag[paramind,cellind]=BGZ5_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ6_lag10_diag[paramind,cellind]=BGZ6_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ7_lag10_diag[paramind,cellind]=BGZ7_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ8_lag10_diag[paramind,cellind]=BGZ8_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ9_lag10_diag[paramind,cellind]=BGZ9_lag10/BGZ_lag10_diag[paramind,cellind]
                    BGZ10_diag[paramind,cellind]=BGZ10/BGZ_lag10_diag[paramind,cellind]
    
                    BGN1_diag[paramind,cellind]=BGN1/BGN_diag[paramind,cellind]
                    BGN2_diag[paramind,cellind]=BGN2/BGN_diag[paramind,cellind]
                    BGN3_diag[paramind,cellind]=BGN3/BGN_diag[paramind,cellind]
                    BGN4_diag[paramind,cellind]=BGN4/BGN_diag[paramind,cellind]
                    BGN5_diag[paramind,cellind]=BGN5/BGN_diag[paramind,cellind]
                    BGN6_diag[paramind,cellind]=BGN6/BGN_diag[paramind,cellind]                
                    BGN7_diag[paramind,cellind]=BGN7/BGN_diag[paramind,cellind]                
    
            
                    ### end cellind loop 
                    ### end paramind loop
            drought_diag=drought1_diag+drought2_diag 
            ###### saving output
            
            for diagvar_name0 in diagvar_list:
                diagvar_name=diagvar_name0+'_diag'
                diagvar=globals()[diagvar_name]
    
            
                ########## save diagnostics to file
                if save_output==1:
                    np.save(exp + '_' + version +'_'+param_name+'_'+diagvar_name+'.npy', diagvar)
    
                ### end diag var
    
        ### end parameter
    ### end variant
### end model
    
    
