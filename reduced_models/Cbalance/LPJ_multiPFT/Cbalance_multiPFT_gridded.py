#!/usr/bin/env python3`bm_inc_ind_LPJ
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 17:05:13 2023

@author: bathiany
"""
import numpy as np
import random
import subprocess
#import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
sys.path.insert(1, '/home/bathiany/Projects/Vegetation_resilience_indicators/reduced_models/Cbalance/functions')
import LPJ_functions
import par
#from ACF import ACF_fun
#from sklearn.linear_model import LinearRegression
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
##np.seterr(invalid='ignore')
#np.seterr(divide='raise', invalid='raise')
np.seterr(divide='ignore', invalid='ignore')
sys.path.insert(1, '/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ')
import assign_properties


#### to do
#    x is it really OK to ignore grasses? it seems so, high corr!
#    x what to do with sudden deaths in LPJ? ignore for now
# x build in intervention in year x, where state variables are modified in some way
# x compute relaxation time with Taylor's code

########## parameters
##exp='sbared0002' # old

# input2 grid
NP=20
NT=20

save_output=1


inicond='LPJ'

versions=("A1N1M1E1A1S0",)  # compare to S1 !

## use all grid cells:
Pindini=0
Tindini=0


exp_driving="sba0056" ## input data from this run
YRS_spinup=500
YRS=YRS_spinup+1000
year_perturb=200

### test
#YRS_spinup=0
##YRS=YRS_spinup+20
#year_perturb=10


# all
#perturb_list=('0010', '0011', '0012', '0013', '0014', '0015', 
#              '0020', '0021', '0022', '0023', '0024', '0025', 
#              '0030', '0031', '0032', '0033', '0034', '0035')
#perturb_list=('0041','0042','0043','0044')


##long stationary run (P0000)
#YRS_spinup=0
#YRS=10000

LPJpath='/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/'

for perturb in perturb_list:    
    
    exp=exp_driving+'_P'+perturb
    
    #######################
    
    YRS_LPJ=YRS+1   # in reality, number of years is YRS. But [0:YRS] actually ends with YRS-1, stupid python
    if (YRS_LPJ>10000):
        YRS_LPJ=10000
    
    if exp_driving == "sba0056" or exp_driving == "sba0042":
        NPFT=8
    else:
        NPFT=1
    
    ### perturbation factors
    # here: affect all PFTs in the same way
    # from 0   
    if perturb == '0000':
        perturbation_factors=[1,1,1,1,1]   # nothing is disturbed (control run)
    elif perturb == '0010':
        perturbation_factors=[0,0,0,0,0]   # all to 0
    elif perturb == '0011':
        perturbation_factors=[0,1,1,1,1]   # all leaf to 0
    elif perturb == '0012':
        perturbation_factors=[1,0,1,1,1]   # all root to 0
    elif perturb == '0013':
        perturbation_factors=[1,1,0,1,1]
    elif perturb == '0014':
        perturbation_factors=[1,1,1,0,1]
    elif perturb == '0015':
        perturbation_factors=[1,1,1,1,0]

    # from below        
    elif perturb == '0020':
        perturbation_factors=[0.5,0.5,0.5,0.5,0.5]
    elif perturb == '0021':
        perturbation_factors=[0.5,1,1,1,1]
    elif perturb == '0022':
        perturbation_factors=[1,0.5,1,1,1]
    elif perturb == '0023':
        perturbation_factors=[1,1,0.5,1,1]
    elif perturb == '0024':
        perturbation_factors=[1,1,1,0.5,1]
    elif perturb == '0025':
        perturbation_factors=[1,1,1,1,0.5]

    # from above    
    elif perturb == '0030':
        perturbation_factors=[2,2,2,2,2]
    elif perturb == '0031':
        perturbation_factors=[2,1,1,1,1]
    elif perturb == '0032':
        perturbation_factors=[1,2,1,1,1]
    elif perturb == '0033':
        perturbation_factors=[1,1,2,1,1]
    elif perturb == '0034':
        perturbation_factors=[1,1,1,2,1]
    elif perturb == '0035':
        perturbation_factors=[1,1,1,1,2]

    # special: exchange PFTs' nind
    elif perturb == '0040':  # expand: PFT 1, 3, 5, 7  # reduce: PFT 2, 4, 6, 8
        perturbation_factors_N=[2.0, 0.5, 2.0, 0.5, 2.0, 0.5, 2.0, 0.5]
    elif perturb == '0041':  # opposite
        perturbation_factors_N=[0.5, 2.0, 0.5, 2.0, 0.5, 2.0, 0.5, 2.0]
    elif perturb == '0042':  # expand: PFT 1, 2, 4, 6, # reduce: PFT 3, 5, 7, 8
        perturbation_factors_N=[2.0, 2.0, 0.5, 2.0, 0.5, 2.0, 0.5, 0.5]   
    elif perturb == '0043':  # opposite
        perturbation_factors_N=[0.5, 0.5, 2.0, 0.5, 2.0, 0.5, 2.0, 2.0]           
    elif perturb == '0044':  # expand: PFT 2, 3, 8,   # reduce: PFT 1, 4, 5, 6, 7
        perturbation_factors_N=[0.5, 2.0, 2.0, 0.5, 0.5, 0.5, 0.5, 2.0]           
    elif perturb == '0045':  # opposite
        perturbation_factors_N=[2.0, 0.5, 0.5, 2.0, 2.0, 2.0, 2.0, 0.5]           
    
    
    
    #### dims are: time, PFT, P, T
    
    ### drivers
    
    file = nc.Dataset(LPJpath+exp_driving+'/pft_nind.nc')
    N_afteryear = file['nind']
    N_afteryear_LPJ=LPJ_functions.unmask(N_afteryear[0:YRS_LPJ,0:NPFT,:,:])
    
    # determine pseudo spatial dimensions:
    NP=np.shape(N_afteryear_LPJ)[2]
    NT=np.shape(N_afteryear_LPJ)[3]
    
    
    file = nc.Dataset(LPJpath+exp_driving+'/pft_nind_beforemort.nc')
    if exp_driving=='sba0042':
        N_beforemort = file['pft_nind_beforemort']  # sba0042
    else:
        N_beforemort = file['pft_nind']   # sba0056, sba0053, sba0044
    N_beforemort_LPJ=LPJ_functions.unmask(N_beforemort[0:YRS_LPJ,0:NPFT,:,:])
    
    
    ## sba0042 # perhaps is wrong? need it before allo!
    if exp_driving=='sba0042' or exp_driving=='sba0044':
        file = nc.Dataset(LPJpath+exp_driving+'/pft_bm_inc_carbon.nc',mode='r')
        bm = file.variables['bm_inc_carbon']
    else:
        file = nc.Dataset(LPJpath+exp_driving+'/pft_bm_inc_carbon_beforeallo.nc',mode='r')
        bm = file.variables['bm_inc_carbon_beforeallo']
    
    bm_inc_LPJ=LPJ_functions.unmask(bm[0:YRS_LPJ,0:NPFT,:,:])
    
    ## compute indiv-based bm_inc while taking care of inactive PFTs.
    #bm_inc_LPJ = np.float32(bm_inc_LPJ)
    #N_beforemort_LPJ = np.float32(N_beforemort_LPJ)
    
    bm_inc_ind_LPJ=bm_inc_LPJ/N_beforemort_LPJ
    bm_inc_ind_LPJ[np.isinf(bm_inc_ind_LPJ)]=0.
    bm_inc_ind_LPJ[np.isnan(bm_inc_ind_LPJ)]=0.
    
    ## Convert masked values to 0 (otherwise get '--')
    #bm_inc_ind_LPJ = np.ma.filled(bm_inc_ind_LPJ, fill_value=0.)
    
    file = nc.Dataset(LPJpath+exp_driving+'/pft_wscal.nc')
    wscal_all = file['pft_wscal']
    wscal_LPJ=LPJ_functions.unmask(wscal_all[0:YRS_LPJ,0:NPFT,:,:])
    
    
    ## time-constant n_woody
    file = nc.Dataset(LPJpath+exp_driving+'/n_woody_timmean.nc')
    n_woody_all = file['pft_n_est']
    n_woody_LPJ=LPJ_functions.unmask(n_woody_all[0,:,:])
    
    
    ### ind pools, LPJ, for inicond
    file = nc.Dataset(LPJpath+exp_driving+'/pft_cleaf.nc')
    L_all = file['C_leaf_pft']
    L_LPJ=LPJ_functions.unmask(L_all[0:YRS_LPJ,0:NPFT,:,:])
    
    file = nc.Dataset(LPJpath+exp_driving+'/pft_croot.nc')
    R_all = file['C_root_pft']
    R_LPJ=LPJ_functions.unmask(R_all[0:YRS_LPJ,0:NPFT,:,:])
    
    file = nc.Dataset(LPJpath+exp_driving+'/pft_csapw.nc')
    S_all = file['C_sapwood_pft']
    S_LPJ=LPJ_functions.unmask(S_all[0:YRS_LPJ,0:NPFT,:,:])
    
    file = nc.Dataset(LPJpath+exp_driving+'/pft_chawo.nc')
    H_all = file['C_heartwood_pft']
    H_LPJ=LPJ_functions.unmask(H_all[0:YRS_LPJ,0:NPFT,:,:])
    
    #file = nc.Dataset(LPJpath+exp_driving+'/pft_cdebt.nc')
    #D_all = file['C_debt_pft']
    #D_LPJ=LPJ_functions.unmask(D_all[0:YRS_LPJ,0:NPFT,:,:])
    
    C_LPJ=L_LPJ+R_LPJ+S_LPJ+H_LPJ
    
    
    file = nc.Dataset(LPJpath+exp_driving+'/pft_height.nc')
    height_all = file['pft_height']
    height_LPJ=LPJ_functions.unmask(height_all[0:YRS_LPJ,0:NPFT,:,:])
    
    #file = nc.Dataset(LPJpath+exp_driving+'/vegc.nc')
    #vegc_all = file['VegC']
    #vegc_LPJ=LPJ_functions.unmask(vegc_all[0:YRS_LPJ,:,:])
    #
    #file = nc.Dataset(LPJpath+exp_driving+'/agb.nc')
    #agb_all = file['AGB']
    #agb_LPJ=LPJ_functions.unmask(agb_all[0:YRS_LPJ,:,:])
    
    
    ###### climate inputs
    #inputpath="/home/bathiany/pik/input2/"
    #file = nc.Dataset(inputpath+'Tair_daily_10000yrs_domain20x20.nc')
    #Prec_all = file['Tair']
    #Prec_LPJ=Prec_all[0:YRS_LPJ,:,:]
    
    ##### time mean inputs
    bm_inc_ind_LPJ_tm=np.mean(bm_inc_ind_LPJ, axis=0)
    
    ## all the same
    #n_woody_LPJ_tm=np.mean(n_woody_LPJ, axis=0)
    #n_woody_LPJ_timmax=np.max(n_woody_LPJ, axis=0)
    wscal_LPJ_tm=np.mean(wscal_LPJ, axis=0)
    
    ## time means only for eval
    #vegc_LPJ_tm=np.mean(vegc_LPJ, axis=0)
    #agb_LPJ_tm=np.mean(agb_LPJ, axis=0)
    
    
    par.set_fixed_parameters_multiPFT()
    
    fL=par.fL
    fR=par.fR
    fS=par.fS
    k_est=par.k_est             # same as est_max, maximum overall sapling establishment rate (indiv/m2)
    fpcmax=par.fpc_tree_max
    
    Nvarfrac=1
    estfrac=1
    est_prescribed=-999
    mort_prescribed=-999
    
    
    def write_data(year_out_index):
        pft_cleaf[year_out_index,:,Pind,Tind]=L
        pft_croot[year_out_index,:,Pind,Tind]=R
        pft_csapw[year_out_index,:,Pind,Tind]=S
        pft_chawo[year_out_index,:,Pind,Tind]=H
        pft_nind[year_out_index,:,Pind,Tind]=N
        fpc[year_out_index,:,Pind,Tind]=FPC
    
    for version in versions:
    
        model=version[0:4]
        adj_switch=int(version[9])
        shuffle=int(version[11])
        #############################
    
        #### OUTPUT
    
        #  full data with time dimension and PFT    
        pft_cleaf=np.zeros((YRS-YRS_spinup,NPFT,NP,NT))
        pft_croot=np.zeros((YRS-YRS_spinup,NPFT,NP,NT))
        pft_csapw=np.zeros((YRS-YRS_spinup,NPFT,NP,NT))
        pft_chawo=np.zeros((YRS-YRS_spinup,NPFT,NP,NT))
        pft_nind=np.zeros((YRS-YRS_spinup,NPFT,NP,NT))
        fpc=np.zeros((YRS-YRS_spinup,NPFT,NP,NT))
    
        # everything else should be computed from this output with ksh scripts
        # so that this model does not have to be run again for more diagnosis
    
        print(exp)
        for Tind in range(Tindini,NT):
            print(Tind, '/', NT)

            for Pind in range(Pindini,NP):
                #print(Pind,Tind)
    
                n_woody=n_woody_LPJ[Pind, Tind]
               
                ###### initialise from full LPJ model in 1979          
                L=L_LPJ[0,:,Pind, Tind]
                R=R_LPJ[0,:,Pind, Tind]
                S=S_LPJ[0,:,Pind, Tind]
                H=H_LPJ[0,:,Pind, Tind]
                #D=D_LPJ[0,:,Pind, Tind]
                D=np.zeros(NPFT)
                N=N_afteryear_LPJ[0,:,Pind, Tind]
                height=height_LPJ[0,:,Pind, Tind]
                FPC=np.zeros(NPFT) + fpcmax # not needed as inicond for the model
                           
                year_out_index=0
                if YRS_spinup==0:
                    write_data(year_out_index)   # updates year_out_index and the output arrays
                    year_out_index+=1
    
                year_out=1979      # 1979 was the initialisation
                for year_index in range(1,YRS):   # year_index counts everything, starting at 1. Last index is YRS-1, but with the initial condition we do have YRS years.
                    
                    if shuffle==0:   # no shuffling in any input
                        year_A=year_index
                        #year_N=year_index
                        #year_M=year_index
                        #year_E=year_index
            
                    elif shuffle==1:  # shuffle all inputs the same way (destroys AC, but keeps corr)
                        #year_A=random.randint(1, YRS_LPJ-1)
                        year_A=random.randint(0, YRS_LPJ)
                        #year_N=year_A
                        #year_M=year_A
                        #year_E=year_A
                        
                    if shuffle==3:    # time const input
                        bm_inc_ind_input = bm_inc_ind_LPJ_tm[:, Pind, Tind]
                        wscal_input=wscal_LPJ_tm[:, Pind, Tind]
                        #n_woody=n_woody_LPJ_timmax[Pind, Tind]
                        #n_woody=n_woody_LPJ[0,Pind, Tind] # see above; only change when time-dependent n_woody is used
    
                    else:   # time-dependent input
                        bm_inc_ind_input = bm_inc_ind_LPJ[year_A, :, Pind, Tind]
                        wscal_input=wscal_LPJ[year_A, :, Pind, Tind]
                                                        
                    ### run the model
                    if model=='A1N1':
                        L, R, S, H, D, N, FPC, height, adj_count, adjN, adjNrel, mort, Nmort, est, estvsN = \
                               LPJ_functions.LPJ_Cbalance_A1N1_multiPFT(L, R, S, H, \
                                  D, N, height, bm_inc_ind_input, \
                                  wscal_input, \
                                  fL, fR, fS, Nvarfrac, k_est, est_prescribed, estfrac, n_woody, \
                                  mort_prescribed, fpcmax, adj_switch, NPFT)
                            
                    else:
                        print("model version not found")
#        
#                    if year_index >= YRS_spinup:
#                        write_data(year_out_index)
#                        year_out_index+=1
#                        #print(L, R, S, H, N, year_out_index)
#                    year_out+=1  # the actual calendar year
                    
                    if year_index == YRS_spinup + year_perturb:
                        ### impose perturbation
                        #L=L*perturbation_factors[0] ## and so on, also works

                        ## special: perturb PFTs in N, each PFT differently
                        if perturb[2]=='4':
                            #print(N)
                            N=N*perturbation_factors_N
                            #print(N)
                        else:
                            ## perturb state
                            state = [L, R, S, H, N]
                            state = [variable * factor for variable, factor in zip(state, perturbation_factors)]
                            L, R, S, H, N = state                            
        

                    if year_index >= YRS_spinup:
                        write_data(year_out_index)
                        year_out_index+=1
                        #print(L, R, S, H, N, year_out_index)
                    year_out+=1  # the actual calendar year
                    

                    
                #### end of time loop
            
        ########## save output to netcdf files
        if save_output==1:
    
            command = ['mkdir -p ' + exp]
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            y_ini= 1979 + YRS_spinup
            y_fin = 1979 + YRS - 1
                
            varlist=('pft_cleaf','pft_croot','pft_csapw','pft_chawo','pft_nind','fpc')
    
            for var in varlist:
                varname=assign_properties.assign_varname(var,'')[0]           
    
                ## copy var from exp_driving
                outfile=exp+'/'+var+'.nc'
                templatefile=LPJpath + exp_driving + '/'+var+'.nc'
                command = ['cdo mulc,0 -sellevel,1/'+str(NPFT)+' -selyear,'+str(y_ini)+'/'+str(y_fin) + ' ' + templatefile + ' ' + outfile]
        
                result = subprocess.run(command, shell=True, capture_output=True, text=True)
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                print("Return Code:", result.returncode)   
                    
                nc_file = nc.Dataset(outfile, mode="a")
        
                data = globals()[var]  # Or use locals() if it's in a local scope
                assert data.shape == (
                    len(nc_file.dimensions["time"]),
                    len(nc_file.dimensions["npft"]),
                    len(nc_file.dimensions["latitude"]),
                    len(nc_file.dimensions["longitude"]),
                ), "Shape mismatch between variable and NetCDF dimensions"
    
                var_nc = nc_file.variables[varname]
                
                # Write data into the variable
                var_nc[:, :, :, :] = data
                
                # Close the NetCDF file
                nc_file.close()
                
        
