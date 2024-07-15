#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:05:54 2022

@author: bathiany
"""

import numpy as np

indices=np.zeros((3,2))

def assign_Tvec(exp):

    expnumber=int(exp[5:7])
    print(exp, expnumber)
    
    if expnumber < 56:  
    #   if exp=='sba0001' or exp=='sba0002' or exp=='sba0003' or exp=='sba0006' or exp=='sba0016' or exp=='sba0022' or exp=='sba0023' or exp=='sba0024' or exp=='sba0025' or exp=='sba0026' or exp=='sba0030' or exp=='sba0031' or exp=='sba0032' or exp=='sba0033' or exp=='sba0034' or exp=='sba0035' or exp=='sba0036' or exp=='sba0042' or exp=='sba0044':
    ## new tropical range, 12-42 °C
        Tvec=[np.round(x,1) for x in np.arange(12,42,2)]
#       indices[0,:]=[10,5]
#       indices[1,:]=[3,5]
#        indices[2,:]=[10,1]

#    elif expnumber...:    # the old tropics, smaller range
#        Tvec=[np.round(x,1) for x in np.arange(17,32,1)]
#        indices[0,:]=[10,10]
#        indices[1,:]=[3,10]
#        indices[2,:]=[10,3]

    elif expnumber >= 56 and expnumber < 80:
        Tvec=[np.round(x,1) for x in np.arange(-8,32,2)]

    elif expnumber >= 80:
        Tvec=[np.round(x,1) for x in np.arange(-999,-998,2)]

    Tvec=np.array(Tvec)        

    return Tvec #, indices


def assign_Pvec(exp):

    expnumber=int(exp[5:7])
    
    if expnumber < 56:  
    #   if exp=='sba0001' or exp=='sba0002' or exp=='sba0003' or exp=='sba0006' or exp=='sba0016' or exp=='sba0022' or exp=='sba0023' or exp=='sba0024' or exp=='sba0025' or exp=='sba0026' or exp=='sba0030' or exp=='sba0031' or exp=='sba0032' or exp=='sba0033' or exp=='sba0034' or exp=='sba0035' or exp=='sba0036' or exp=='sba0042' or exp=='sba0044':
    ## always use new range, 12-42 °C
        Pvec=[np.round(x,1) for x in np.arange(0,3000,200)]
    elif expnumber >= 56 and expnumber < 80:
        Pvec=[np.round(x,1) for x in np.arange(100,2100,100)]        
    elif expnumber > 80:
        Pvec=[np.round(x,1) for x in np.arange(0,4000,200)]        
        
    Pvec=np.array(Pvec)
    return Pvec


def assign_PSvec(exp):

    expnumber=int(exp[5:7])
    
    if expnumber <= 80:
        PSvec=[np.round(x,1) for x in np.arange(-999,-998,2)]
    elif expnumber > 80:
        PSvec=[np.round(x,1) for x in np.arange(0,4950,4950/20)]
    PSvec=np.array(PSvec)
    return PSvec

def assign_PSfactvec(exp):

    expnumber=int(exp[5:7])
    
    if expnumber <= 80:
        PSfactvec=[np.round(x,1) for x in np.arange(-999,-998,2)]
    elif expnumber > 80:
        PSfactvec=[np.round(x,1) for x in np.arange(0,3,3/20)]

    PSfactvec=np.array(PSfactvec)
    return PSfactvec


def assign_varname(var_case,type):

    #if varbase.find('_') != -1:
    #    var_case=varbase.split('_')[0]
    #else:
    #    var_case=varbase

    #print(var_case)

    cmap_var="viridis"
    
    unit=''
    minval=-0.3
    maxval=1
    if var_case=='mnpp':
        varname='NPP'
        unit='[gC/m2/mon]'
        
    elif var_case=='npp':
        varname='NPP'
        unit='[gC/m2/yr]'
        minval=0
        maxval=1500

        if type=='timstd':
            minval=0
            maxval=150

    elif var_case=='mgpp':
        varname='GPP'
        unit='[gC/m2/mon]'
        minval=0
        maxval=100

    elif var_case=='gpp':
        varname='GPP'
        unit='[gC/m2/yr]'        
        minval=0
        maxval=2000
        
    elif var_case=='mra':
        varname='GPP'
        minval=0
        maxval=200
        
    elif var_case=='vegc' or var_case=='Dvegc':
        varname='VegC'
        minval=0
        maxval=30000
        unit="gC/m²"

        if type=='timstd' or type == "SDy":
            minval=0
            maxval=600

#        elif type[0:8]=='recovery':
#            minval=0
#            maxval=500
            
        
    elif var_case=='pft_vegc' or var_case=='Dpft_vegc':
        varname='vegC_pft'
        minval=0
        maxval=28000
        unit="gC/m²"

        if type=='timstd':
            minval=0
            maxval=1000

    elif var_case=='pft_P':
        varname='vegC_pft'
        minval=0
        maxval=1
        unit="ind/m2"
        
    elif var_case=='vegc2' or var_case=='Dvegc2':
        varname='C_leaf_pft'
        minval=0
        maxval=1.5e6

    elif var_case=='agb' or var_case=='Dagb':
        varname='AGB'
        unit="gC/m²"
        minval=0
        maxval=24000
        
    elif var_case=='fpc' or var_case=='fpc_TropBL' or var_case=='fpc_woody':
        varname='FPC'
        unit="m²/m²"
        minval=0
        maxval=1
        
    elif var_case=='pft_mort' or var_case=='pft_mort_diag':
        varname='mortality'
        minval=0
        maxval=0.03
        unit="1/yr"
        
    elif var_case=='pft_cleaf':
        varname='C_leaf_pft'
        unit="gC/m²"
        minval=0
        maxval=20000
        unit="gC/ind"

    elif var_case=='cleaf' or var_case=='Dcleaf':
        varname='C_leaf_pft'
        unit="gC/m²"
        minval=0
        maxval=30000
        unit="gC/ind"        

    elif var_case=='pft_csapw' or var_case=='csapw' or var_case=='Dsapw':
        varname='C_sapwood_pft'
        unit="gC/m²"
        minval=0
        maxval=500000
        unit="gC/ind"        

    elif var_case=='croot' or var_case=='Dcroot':
        varname='C_root_pft'
        unit="gC/m²"
        minval=0
        maxval=30000
        unit="gC/ind"

    elif var_case=='pft_croot':   
        varname='C_root_pft'
        unit='[gC/m2]'
        minval=0
        maxval=20000
    
    elif var_case=='pft_chawo' or var_case=='chawo' or var_case=='Dchawo':
        varname='C_heartwood_pft'
        minval=0
        maxval=500000
        unit="gC/ind"
    
    elif var_case=='pft_cleaf_frac':
        varname='C_leaf_pft'
        minval=0
        maxval=0.125 #0.25
        unit="gC/gC"        
        
    elif var_case=='pft_croot_frac':
        varname='C_root_pft'
        minval=0
        maxval=0.125 #0.25
        unit="gC/gC" 

    elif var_case=='pft_csapw_frac':
        varname='C_sapwood_pft'
        minval=0.1 #0.2
        maxval=0.5 #
        unit="gC/gC" 

    elif var_case=='pft_chawo_frac':
        varname='C_heartwood_pft'
        minval=0.4 #0.5
        maxval=0.75 #0.8
        unit="gC/gC"

    elif var_case=='pft_SH_frac':
        varname='C_sapwood_pft'
        unit="gC/gC"
        minval=0.8
        maxval=0.95


    elif var_case=='litc' or var_case=='Dlitc':
        varname='LitC'
        unit='[gC/m2]'
        minval=0
        maxval=10000
        unit="gC/m2"
        
    elif var_case=='cliving' or var_case=='Dcliving':
        varname='C_heartwood_pft'
        unit="gC/m²"
        minval=0
        maxval=24000
        
    elif var_case=='pft_chawoPLUSpft_sapw':
        unit="gC/m²"
        varname='C_heartwood_pft'
        minval=0
        maxval=300000

    elif var_case=='pft_chawoperm2PLUSpft_sapwperm2':
        unit="gC/m²"
        varname='vegC_pft'
        minval=0
        maxval=600000
        
    elif var_case=='chawoperm2PLUSlit':
        unit="gC/m²"
        varname='vegC_pft'
        minval=0
        maxval=35000
        
    elif var_case=='prod_turnover':
        varname='prod_turnover'
        
    elif var_case=='flux_estab':
        varname='estabc'
        
    elif var_case=='agb_tree':
        varname='agb_tree'
        unit="gC/m²"
        
    elif var_case=='cleafperm2' or var_case=='pft_cleafperm2' or var_case=='crootperm2' or var_case=='pft_crootperm2' or var_case=='csapwperm2' or var_case=='pft_csapwperm2' or var_case=='chawoperm2' or var_case=='pft_chawoperm2':
        unit="gC/m²"
        varname='nind'
        minval=0
        maxval=35000
        
    elif var_case=='pft_vegcperind':
        unit='[gC/ind]'
        varname='vegC_pft'
        minval=0
        maxval=1000000
        
    elif var_case=='pft_vegcfrac':
        unit='[gC/ind]'
        varname='vegC_pft'
        minval=0
        maxval=1

    elif var_case=='pft_frac_csapw':
        unit=''
        varname='nind'
        minval=0.2
        maxval=0.6
        
    elif var_case=='pft_frac_chawo':
        unit=''
        varname='nind'
        minval=0.5
        maxval=0.8
                
#    elif var_case=='pft_est':
#        unit='ind/m2'
#        varname='est'
#        minval=0
#        maxval=0.01
        
    elif var_case=='bm_inc_carbon':
        unit="gC/m²"
        varname="bm_inc_carbon"
        minval=0
        maxval=0
    
    elif var_case=='pft_nind' or var_case=='nind':
        unit='ind/(m²yr)'
        varname='nind'
        minval=0
        maxval=0.09
        
#    elif var_case=='pft_est_diag':
#        unit='ind/m2'
#        varname='nind'
#        minval=0
#        maxval=2
        
    elif var_case=='pft_tinc_ind_cleaf':
        unit="gC/m²"
        varname='tinc_ind_cleaf'
        minval=-2
        maxval=300
        
    elif var_case=='pft_tinc_ind_croot':
        unit="gC/m²"
        varname='tinc_ind_croot'
        minval=-2
        maxval=300
        
    elif var_case=='pft_tinc_ind_csapw':
        unit="gC/m²"
        varname='tinc_ind_csapw'
        minval=-2
        maxval=300
                
    elif var_case=='pft_tinc_ind_chawo':
        unit="gC/m²"
        varname='tinc_ind_chawo'
        minval=-2
        maxval=300
    elif var_case=='pft_n_est':
        varname='pft_n_est'
        minval=0
        maxval=4
        #cmap_var="rainbow"        
        
    elif var_case=='Shannon_fpc':
        varname='FPC'
        unit=''
        minval=0
        maxval=1.3
    elif var_case=='Shannon_pft_nind':
        varname='pft_nind'
        unit=''
        minval=0
        maxval=1.3
        
    elif var_case=='pft_wscal':
        varname='pft_wscal'
        unit=''
        minval=0.0
        maxval=1.0

    elif var_case=='pft_bm_inc_ind':
        varname='bm_inc_carbon_beforeallo'
        unit='gC/ind'
        minval=0.0
        maxval=1.0
    
    elif var_case=='vegcperind':
        varname='VegC'
        unit=''
        minval=0.0
        maxval=1.0

    elif var_case=='est' or var_case=='pft_est':
        varname='establishment'
        unit='ind/(m2yr)'
        minval=0.0
        maxval=0.12

    elif var_case=='adjN':
        varname='pft_nind'
        unit=''
        minval=0.0
        maxval=0.0016

    elif var_case=='adjNrel':
        varname='pft_nind'
        unit=''
        minval=0.0
        maxval=0.018

    elif var_case=='height' or var_case=='pft_height':
        varname='pft_height'
        unit=''
        minval=0.0
        maxval=25  


    else:
        print('variable not found')

    return varname, unit, minval, maxval, cmap_var