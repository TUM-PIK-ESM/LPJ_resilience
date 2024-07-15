#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 17:05:13 2023

@author: bathiany
"""
import numpy as np
import os
import sys
sys.path.insert(1, '/home/bathiany/Projects/Amazon_EWS/veg_EWS_reducedmodels/Cbalance/functions')
import par
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
np.seterr(invalid='ignore')

########## parameters
dpi=500
exp="sba0084"      # exp from which data is read in
PFT=1

Nparams=7

Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]
Pvec=np.array(Pvec)                 # in mm/year

#Nprec=19   # 
cellind_ini=1   # remove P=0, starting from 1 (200 mm/yr, i.e. cellind=1)
cellind_fin=13
Pvec=Pvec[cellind_ini:cellind_fin] 

Nprec=cellind_fin-cellind_ini

############ paper

#### show that varfrac even reverses pattern when using fixed N
#models=("A1N3",)
#variants=("M0E1A1",)
#shuffle=1
#parameters=("varfrac",)
#diagvar_list=("AC10_Cperm2","AC30_Cperm2", "drought","SHfrac")


## for all versions plot
#versions=("A1N1M1E1A1S1","A1N1M0E1A1S1","A1N3M0E1A1S1","A2N1M0E1A1S1")
#models=("A1N1","A1N3","A2N1")
models=("A1N1",)
variants=("M0E1A1",)
shuffle=1
parameters=("Nvarfrac",)
#parameters=("varfrac",)
diagvar_list=("AC10_Cperm2","AC30_Cperm2")


##### separation of terms
#diagvar_list=("BGZ1_lag10","BGZ2_lag10",\
#              "BGZ3_lag10","BGZ4_lag10",\
#              "BGZ5_lag10","BGZ6_lag10",\
#              "BGZ7_lag10","BGZ8_lag10",\
#              "BGZ9_lag10","BGZ10",\
#              "BGN1","BGN2","BGN3","BGN4","BGN5","BGN6","BGN7")


### role of Ndyn
#models=("A2N1","A1N1",)
#models=("A1N1",)
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
#parameters=("varfrac",)
#diagvar_list=("AC10_Cperm2",)

## ACheatmaps, variance scaling
## Allocation and variance: why does AC decrease with less variance (at dry cells in N1 and all cells in N0)?
models=("A1N3",)
variants=("M0E1A1",)
shuffle=1
parameters=("varfrac",)
diagvar_list=("AC10_Cperm2","drought","SHfrac")

##diagvar_list=("AC10_Cperm2","height","Hfrac","AHfrac","SHfrac","ASHfrac",)
#diagvar_list=("AC10_Cperm2",)
##diagvar_list=("AC10_Cperm2","adjNrel","drought1","drought2","drought","mort","SHfrac","Hfrac","VarCperm2")
##diagvar_list=("AC10_Cperm2","adjNrel","drought","mort","SHfrac","Hfrac")
##diagvar_list=("AC10_Cperm2","drought","SHfrac","est", "corr_Nest")
##diagvar_list=("corr_Nest",)







### A0: prescribe allocated fractions:

#
#models=("A0N0",)
##models=("A0N3","A1N3")
#variants=("M0E1A1",)
#shuffle=1
#parameters=("varfrac",)
##diagvar_list=("AC10_Cperm2","SHfrac","height")
#diagvar_list=("AC10_Cperm2",)
#
#
##models=("A0N0","A0N1","A0N3")
##variants=("M1E1A1","M0E1A1",)
##shuffle=1
##parameters=("varfrac","ASHfrac")
##diagvar_list=("AC10_Cperm2","Cperm2","SHfrac")
#models=("A0N3",)
#variants=("M0E1A1",)
#shuffle=1
#parameters=("varfrac",)
#diagvar_list=("AC10_Cperm2","Cperm2","SHfrac",)


#models=("A0N3",)
#variants=("M0E1A1",)
#shuffle=1
##parameters=("varfrac",)
##parameters=("AHfrac","ASHfrac")   # both has no effect!
#parameters=("ASHfrac",)
##diagvar_list=("AC10_Cperm2","Hfrac","SHfrac","height")
#diagvar_list=("AC10_Cperm2",)

# Quark
#models=("A0N3",)
#variants=("M0E1A1",)
#shuffle=1
#parameters=("ASHfrac",)
#diagvar_list=("AC10_Cperm2","SHfrac")




## understand variance of biomass: sensitivity to prec variability versus other effect:
#models=("A1N1",)
#variants=("M0E1A1","M1E1A1")
#shuffle=1
#parameters=("varfrac",)
#diagvar_list=("AC10_Cperm2",)




############ parameters to be varied
par.set_fixed_parameters(PFT)


###### plots and saving output


for model in models:
    for variant in variants:
        version=model+variant+'S'+str(shuffle)
        adj_switch=int(variant[5])

        
        for parameter in parameters:
        
            ## set parameter values
            par.set_variable_parameters(Nparams, parameter, paramind=0)
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
            
            for diagvar_name in diagvar_list:
    
                filename=exp + '_' + version + '_'+parameter+'_vec_'+diagvar_name+'_diag.npy'
                if os.path.isfile(filename):
                    data=np.load(filename)
 
#            ### masking
#                    asym_diag=np.load(exp + '_' + version + '_'+parameter+'_vec_'+'asym_diag.npy')
#                    ACFmin_diag=np.load(exp + '_' + version + '_'+parameter+'_vec_'+'ACFmin_diag.npy')
#                    
#                    for paramind in range(0,Nparams):
#                        for cellind in range(0,Nprec):
#                            if asym_diag[paramind,cellind]>4 or ACFmin_diag[paramind,cellind]<-0.6:
#                            # if (np.min(wscal_input[:, cellind+1]) < 0 or np.min(bm_inc_ind_input[:, cellind+1]) < 0  or np.min(N_afteryear_input[:, cellind+1])<0):
#                            ######## Due to rescaling, some inputs can have become < 0
#                            # There, need to mask the diagnosed output
#                            #if (version[0:2]=="A0") and ( np.min(AL_input[year_A, cellind+1] < 0 or np.min(AR_input[year_A, cellind+1] < 0)) :
#                        #                        if (NstdperN_diag[paramind,cellind]>0.05 or Cstd_diag[paramind,cellind]>10000 or Cperm2std_diag[paramind,cellind]>500):    
##                    if asym_diag[paramind,cellind]>4 or ACFmin_diag[paramind,cellind]<-0.6:
#                                data[paramind,cellind]='NaN'


                    plt.figure()

                    var1_name_extended=parameter+'_vec'
                    var1=globals()[var1_name_extended]               # ts is the axis (parameter)
                    var2=data[:,cellind_ini:cellind_fin]                           # this is the color (AC or diagnostic)
                    fig, ax = plt.subplots(figsize=(10,10))

                    if diagvar_name[0:2] == "AC":
#                        vmin=-0.3
#                        vmax=0
                        vmin=0.
                        vmax=1.0
                        #if parameter=='Nvarfrac':
                        #    vmin=0.
                        #    vmax=1.0
                        #elif parameter=='varfrac' and 'model' == 'A1N1':
                        #    vmin=0
                        #    vmax=0.6                  
                        
                        
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
                    elif diagvar_name == "N":
                        vmin=0.1
                        vmax=5
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
                    elif diagvar_name == "adjNrel":
                        vmin=0
                        vmax=0.01
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
                    elif diagvar_name == "estvsN":
                        vmin=0.0152
                        vmax=0.0153
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax) 
                    elif diagvar_name == "AHfrac":
                        vmin=0
                        vmax=0.35
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
                    elif diagvar_name == "ASHfrac":
                        vmin=0
                        vmax=0.35
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
                    elif diagvar_name == "Nstd":
                        vmin=0.0
                        vmax=0.1
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)                            
                    elif diagvar_name == "Cstd":
                        vmin=0.0
                        vmax=3000
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
                    elif diagvar_name == "C":
                        vmin=0.0
                        vmax=500000
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
                    elif diagvar_name == "Cperm2std":
                        vmin=0.0
                        vmax=600
                        Panel=plt.imshow(var2, cmap="viridis", vmin=vmin, vmax=vmax)
    
                    elif diagvar_name[0:3] == "BGZ":
    #                    vmin=-75000
    #                    vmax=75000
                        if diagvar_name == "BGZ_lag10":
                            vmin=-7500
                            vmax=7500
                        else: 
                            vmin=-10
                            vmax=10                 
                        Panel=plt.imshow(var2, cmap="coolwarm", vmin=vmin, vmax=vmax)
    
                    elif diagvar_name[0:3] == "BGN":
                        if diagvar_name == "BGN":
                            vmin=-100000
                            vmax=100000
                        else:
                            vmin=-10
                            vmax=10
                        Panel=plt.imshow(var2, cmap="coolwarm", vmin=vmin, vmax=vmax)                    
    
                    elif diagvar_name[0:4] == "corr":
                        vmin=-1
                        vmax=1
                        Panel=plt.imshow(var2, cmap="coolwarm", vmin=vmin, vmax=vmax)                    
                        
                    else:
                        Panel=plt.imshow(var2, cmap="viridis")

                    fontsize_ticks=24
                    fontsize_title=20
                    fontsize_labels=30
                    ax.invert_yaxis()


                    ax.set_yticks(np.arange(0,Nparams))
                    ax.set_xticks(np.arange(0,Nprec))
                    ax.set_yticklabels(var1, fontsize=fontsize_ticks, rotation=0)
                    ax.set_xticklabels(Pvec, fontsize=fontsize_ticks, rotation=45)
                    
                    cbar = plt.colorbar(Panel, fraction=0.03)


                    if diagvar_name=='AC10_Cperm2':
                        cbar.set_label('autocorrelation (10yr)', rotation=90, fontsize=fontsize_labels)
                    elif diagvar_name=='AC30_Cperm2':
                        cbar.set_label('autocorrelation (30yr)', rotation=90, fontsize=fontsize_labels)
                    elif diagvar_name=='drought':
                        cbar.set_label('drought occurrence', rotation=90, fontsize=fontsize_labels)
                    elif diagvar_name=='SHfrac':
                        cbar.set_label('$(S+H)/C_{ind}$', rotation=90, fontsize=fontsize_labels)
                    else:
                        cbar.set_label(diagvar_name, rotation=90, fontsize=fontsize_labels)
                                
                    
                    cbar.ax.tick_params(labelsize=fontsize_ticks)
                    if parameter == 'Nvarfrac':
                        ax.set_ylabel('$f_N$', fontsize=fontsize_labels)
                    elif parameter == "varfrac":
                        ax.set_ylabel('$f_I$', fontsize=fontsize_labels)    
                    else:
                        ax.set_ylabel(parameter, fontsize=fontsize_labels)
                                        
                    ax.set_xlabel("Precipitation [mm/yr]", fontsize=fontsize_labels)
                    file_fig=exp+'_' + version + '_Prec_vs_'+parameter+'_vs_'+diagvar_name+'_heatmap.png'
                    #plt.title(exp+'_' + version, fontsize=fontsize_title)
                    plt.savefig(file_fig, dpi=dpi, bbox_inches='tight')
                    plt.show()
                else:
                    print('File does not exist: ', filename)