#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 18:49:41 2023

@author: bathiany
"""

import numpy as np
import par

def set_fixed_parameters(PFT):
    global allom1, allom2, allom3 
    global fL_default, fR_default, fS_default
    global vscal,  lmro_ratio, lmro_offset
    global reinickerp, k_latosa, wooddens
    global k_mort, mort_max, height_max, fpc_tree_max_default
    global MAX_ITER, epsilon, CDEBT_MAXLOAN_DEFICIT, CDEBT_MAXLOAN_MASS, NSEG
    global EPSILON_BISECT_x, EPSILON_BISECT_y, MAXITER_BISECT
    global lightextcoeff, sla, CAmax
    global Lsapl, Rsapl, Ssapl, Hsapl, Csapl
    global k_est_default
    
    ## universal model params
    allom1=100.0     # no unit
    allom2=40.0      # no unit
    allom3=0.67      # /*0.5*/   no unit
    reinickerp=1.6   # no unit
    k_latosa=6e3     # no unit  , from par/pft.js
    wooddens=2e5     # gC/m3    , from include/tree.h
    k_mort=0.2       # no unit
    mort_max=0.03       # 1/yr max mortality rate
    height_max=100.0    # maximum height of trees in m
    
    fpc_tree_max_default=0.95   # m2/m2  if set to 0.9, everything already goes crazy
    
    MAX_ITER=10
    epsilon=1.0E-7
    
    lmro_ratio=1.0    ## this is the same for all trees! leaf mass to root mass, from par/pft.js
    lmro_offset=0.5   # from par/pft.js
    
    CDEBT_MAXLOAN_DEFICIT=0.8 # maximum loan as a fraction of deficit  # in allocation_tree.c
    CDEBT_MAXLOAN_MASS=0.2    # maximum loan as a fraction of (sapwood-cdebt)  # in allocation_tree.c
    NSEG=20 # number of segments (parameter in numerical methods)
    #CDEBT_PAYBACK_RATE=0.2     # in turnover_tree.c
    
    ### parameters for leftmostzero
    ## default from LPJ:
    EPSILON_BISECT_x=0.001 #accuracy in x
    EPSILON_BISECT_y=1.0e-10 #accuracy in y           
    MAXITER_BISECT=40 # max iter

    vscal=1  # no nitrogen limitation
    
    k_est_default=0.12
    
    if PFT==1:
        ## TrBE
        fL_default=0.5          # 1/2 in 1/yr
        fR_default=0.5          # 1/2 in 1/yr
        fS_default=0.03333333   # 1/30 in 1/yr
        lai_sapl=1.5    # m2/m2
        wood_sapl=1.2   #
        #longevity=1.6     
        CAmax=25.0      # m2
        lightextcoeff=0.5
        sla=0.01986     # m2/gC => 500g/m2   in m2/m2


    Lsapl=np.power(lai_sapl*allom1*np.power(wood_sapl,reinickerp)*np.power(4.0*sla/np.pi/k_latosa,reinickerp*0.5)/sla,2.0/(2.0-reinickerp));

    stemdiam=wood_sapl*np.sqrt(4.0*Lsapl*sla/np.pi/k_latosa)
    height_sapl=allom2*np.power(stemdiam,allom3)

    Ssapl=wooddens*height_sapl*Lsapl*sla/k_latosa
    Hsapl=(wood_sapl-1.0)*Ssapl
    Rsapl=(1.0/lmro_ratio)*Lsapl
    Csapl=Lsapl+Rsapl+Ssapl+Hsapl


def set_fixed_parameters_multiPFT():
    global allom1, allom2, allom3 
    #global fL_default, fR_default, fS_default
    global fL, fR, fS
    global vscal,  lmro_ratio, lmro_offset
    global reinickerp, k_latosa, wooddens
    global k_mort, mort_max, height_max, fpc_tree_max
    global MAX_ITER, epsilon, CDEBT_MAXLOAN_DEFICIT, CDEBT_MAXLOAN_MASS, NSEG
    global EPSILON_BISECT_x, EPSILON_BISECT_y, MAXITER_BISECT
    global lightextcoeff, sla, CAmax
    global Lsapl, Rsapl, Ssapl, Hsapl, Csapl
    global k_est
    
    ## universal model params
    allom1=100.0     # no unit
    allom2=40.0      # no unit
    allom3=0.67      # /*0.5*/   no unit
    reinickerp=1.6   # no unit   
    k_latosa=6e3     # no unit  , from par/pft.js
    wooddens=2e5     # gC/m3    , from include/tree.h
    k_mort=0.2       # no unit
    mort_max=0.03       # 1/yr max mortality rate
    height_max=100.0    # maximum height of trees in m
    
    fpc_tree_max=0.95   # m2/m2  if set to 0.9, everything already goes crazy
    
    MAX_ITER=10
    epsilon=1.0E-7
    
    ## these two numbers are the same for all 8 tree types! from par/pft.js
    lmro_ratio=1.0    ## lmro=leaf mass to root mass
    lmro_offset=0.5   
    
    CDEBT_MAXLOAN_DEFICIT=0.8 # maximum loan as a fraction of deficit  # in allocation_tree.c 
    CDEBT_MAXLOAN_MASS=0.2   # maximum loan as a fraction of (sapwood-cdebt)  # in allocation_tree.c 
    NSEG=20 # number of segments (parameter in numerical methods)
    #CDEBT_PAYBACK_RATE=0.2     # in turnover_tree.c 
    
    ### parameters for leftmostzero
    ## default from LPJ:
    EPSILON_BISECT_x=0.001 #accuracy in x
    EPSILON_BISECT_y=1.0e-10 #accuracy in y           
    MAXITER_BISECT=40 # max iter

    vscal=1  # no nitrogen limitation
    
    k_est=0.12

    ### The 8 tree types:    
    # tropical broadleaved evergreen tree
    # tropical broadleaved raingreen tree
    # temperate needleleaved evergreen tree
    # temperate broadleaved evergreen tree
    # temperate broadleaved summergreen tree
    # boreal needleleaved evergreen tree
    # boreal broadleaved summergreen tree
    # boreal needleleaved summergreen tree    
    
##   # PFT-specific params for all 8 trees (no grasses)
#    fL_default=[0.5, 1, 0.25, 1, 1, 0.25, 1, 1]
#    fR_default=[0.5, 1, 0.25, 1, 1, 0.25, 1, 1]
#    fS_default=[0.0333333333, 0.033333333333, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
    fL=[0.5, 1, 0.25, 1, 1, 0.25, 1, 1]
    fR=[0.5, 1, 0.25, 1, 1, 0.25, 1, 1]
    fS=[0.0333333333, 0.033333333333, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
    lai_sapl=[1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]
    wood_sapl=[1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
    CAmax=[25, 15, 15, 15, 25, 20, 20, 20]
    lightextcoeff=[0.5, 0.5, 0.4, 0.5, 0.6, 0.4, 0.5, 0.6]
    sla=[0.01986, 0.03233, 0.01049, 0.01986, 0.03233, 0.01049, 0.03233, 0.02118]

## PFT no 8:
    #"aphen_min" : 10,
    #"aphen_max" : 200,
# rest:
    #define APHEN_MAX 245
    #define APHEN_MIN 60    /* minimum aphen for cold-induced senescence */

    NPFT_sapl=8
    Lsapl=np.zeros(NPFT_sapl)
    Rsapl=np.zeros(NPFT_sapl)
    Ssapl=np.zeros(NPFT_sapl)
    Hsapl=np.zeros(NPFT_sapl)
    Csapl=np.zeros(NPFT_sapl)
    
    for PFTind in range(0,8):

 ##     Lsapl=np.power(lai_sapl*allom1*np.power(wood_sapl,reinickerp)*np.power(4.0*sla/np.pi/k_latosa,reinickerp*0.5)/sla,2.0/(2.0-reinickerp));
        Lsapl[PFTind]=np.power(lai_sapl[PFTind]*allom1*np.power(wood_sapl[PFTind],reinickerp)*np.power(4.0*sla[PFTind]/np.pi/k_latosa,reinickerp*0.5)/sla[PFTind],2.0/(2.0-reinickerp));
    
        stemdiam=wood_sapl[PFTind]*np.sqrt(4.0*Lsapl[PFTind]*sla[PFTind]/np.pi/k_latosa)
        height_sapl=allom2*np.power(stemdiam,allom3)
    
        Ssapl[PFTind]=wooddens*height_sapl*Lsapl[PFTind]*sla[PFTind]/k_latosa
        Hsapl[PFTind]=(wood_sapl[PFTind]-1.0)*Ssapl[PFTind]
        Rsapl[PFTind]=(1.0/lmro_ratio)*Lsapl[PFTind]

    Csapl=Lsapl+Rsapl+Ssapl+Hsapl




def set_variable_parameters(Nparams, parameter, paramind):
    global Nvarfrac_vec, varfrac_vec, fL_vec, fR_vec, fS_vec, mort_vec, n_woody_vec
    global AHfrac_vec, ASHfrac_vec, fpcmax_vec, est_vec, k_est_vec, estfrac_vec
    global Nmean_vec, Nstd_vec, NAC_vec
    global Nvarfrac, varfrac, fL, fR, fS, mort, n_woody
    global AHfrac, ASHfrac, fpcmax, est, k_est, estfrac
    global mort_prescribed_param, est_prescribed_param, param_name
    global Nmean, Nstd, NAC

    ############ parameters to be varied
    if  Nparams==7:
        
        Nvarfrac_vec=[0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5]   # default: 1
        
        varfrac_vec=[0.01, 0.05, 0.1, 0.3, 0.6, 1, 1.5]    # orig range, leads to problems
        #varfrac_vec=[0.2,0.4,0.6,0.8,1,1.2,1.4]    # medium reduced range
        #varfrac_vec=[0.4,0.5,0.6,0.7,0.8,0.9,1,]    # very reduced range
        
        fL_vec=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]   # default:0.5
        fR_vec=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]   # default:0.5
        fS_vec=[0.01,0.05,0.1,0.15,0.2,0.3,0.4]     # default:0.0333333
        mort_vec=[0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017]
        AHfrac_vec=[0.05,0.1,0.15,0.2,0.25,0.3,0.35]   # fixed fract for allocation
        ASHfrac_vec=[0.05,0.1,0.15,0.2,0.25,0.3,0.35]   # fixed fract for allocation
        #fpc2_vec=[0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]
        #fpcmax_vec=[0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0 ]   # default: 0.95
        #fpcmax_vec=[0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0]   # default: 0.95
  
    
        fpcmax_vec=[0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]   # default: 0.95

    
    
        n_woody_vec=[1, 2, 3, 4, 5, 6, 7]       # default: 1
        k_est_vec=[0.012,0.03,0.075,0.12,0.2,0.5,1.0]   # default: 0.12
        est_vec=[0.0015,0.01,0.02,0.035,0.05,0.065,0.08]   
        estfrac_vec=[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]   # default: 1

        Nmean_vec=[0.06, 0.1, 0.25, 0.5, 1, 2, 3]
        Nstd_vec=[0, 0.015, 0.03, 0.045, 0.06, 0.075, 0.09]
        NAC_vec=[0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9]


    elif Nparams==1:
        
        Nvarfrac_vec=[1]  
        
    else:
        print("unexpected number of parameter values - need to extend list")


    # set all of the paramaters that are not varied to default value

    # parameters already present in original LPJ
    fL=par.fL_default
    fR=par.fR_default
    fS=par.fS_default
    k_est=par.k_est_default     # same as est_max, maximum overall sapling establishment rate (indiv/m2)
    fpcmax=par.fpc_tree_max_default     # same as est_max, maximum overall sapling establishment rate (indiv/m2)

    ## default parameters, newly introduced to the model
    Nvarfrac=1
    varfrac=1
    mort_prescribed_param=-999
    est_prescribed_param=-999
    estfrac=1

    #if parameter=="fpc2":
    #    n_woody=2
    #else:
    #    n_woody=1
    n_woody=1

    AHfrac=-999
    ASHfrac=-999


    ### here set parameter
    if parameter=="Nvarfrac":
        param_name="Nvarfrac_vec"
        Nvarfrac=Nvarfrac_vec[paramind]
    elif parameter=="varfrac":
        param_name="varfrac_vec"
        varfrac=varfrac_vec[paramind]
    elif parameter=="mort":
        param_name="mort_vec"
        mort_prescribed_param=mort_vec[paramind]
    elif parameter=="fS":
        param_name="fS_vec"
        fS=fS_vec[paramind]
    elif parameter=="fL":
        param_name="fL_vec"
        fL=fL_vec[paramind]
    elif parameter=="fR":
        param_name="fR_vec"
        fR=fR_vec[paramind]
    elif parameter=="n_woody":
        param_name="n_woody_vec"
        n_woody=n_woody_vec[paramind]
    elif parameter=="k_est":
        param_name="k_est_vec"
        k_est=k_est_vec[paramind]
    elif parameter=="AHfrac":
        param_name="AHfrac_vec"
        AHfrac=AHfrac_vec[paramind]
    elif parameter=="ASHfrac":
        param_name="ASHfrac_vec"
        ASHfrac=ASHfrac_vec[paramind]
#                    elif parameter=="fpc2":
#                        param_name="fpc2_vec"
#                        fpc2=fpc2_vec[paramind]
    elif parameter=="fpcmax":
        param_name="fpcmax_vec"
        fpcmax=par.fpcmax_vec[paramind]
    elif parameter=="est":
        param_name="est_vec"
        est_prescribed_param=est_vec[paramind]
    elif parameter=="estfrac":
        param_name="estfrac_vec"
        estfrac=estfrac_vec[paramind]
    elif parameter=="Nmean":
        param_name="Nmean_vec"
        Nmean=Nmean_vec[paramind]
    elif parameter=="Nstd":
        param_name="Nstd_vec"
        Nstd=Nstd_vec[paramind]
    elif parameter=="NAC":
        param_name="NAC_vec"
        NAC=NAC_vec[paramind]
    else:
        print("no parameter chosen")