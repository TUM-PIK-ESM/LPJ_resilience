#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 18:36:57 2023

@author: bathiany
"""
import numpy as np
import par


### read LPJ output files are masked arrays,
## convert to normal ones:
def unmask(data):
    data2 = np.float32(data)
    data_unmasked = np.ma.filled(data2, fill_value=0.)
    return data_unmasked
    



###### individual modules (processes)


### missing:
# light competition

# phen computed in lpj/phenology_gsi.c, 
# is used in lpj/water_stressed.c, tree/turnover_daily_tree.c, tree/npp_tree.c, 
# tree/lai_tree.c (see below), lpj/daily_natural.c
# landuse/daily_setaside.c
# landuse/daily_biomass_tree.c
#  and others

## function fpc_tree

## old (single PFT)
#def fpc_tree(L,N,CAind):
#    if CAind>0.0:
#        LAIind=lai_tree(L, CAind)
#        FPCind=1.0-np.exp(-par.lightextcoeff*LAIind)
#        FPC_PFT=CAind*N*FPCind
#    else:
#        FPC_PFT=0.0
#    return FPC_PFT

def fpc_tree(L,N,CAind, sla, lightextcoeff):
    if CAind>0.0:
        LAIind=lai_tree(L, CAind, sla)
        FPCind=1.0-np.exp(-lightextcoeff*LAIind)
        FPC_PFT=CAind*N*FPCind
    else:
        FPC_PFT=0.0
    return FPC_PFT


## old (single PFT)
#def lai_tree(L, CAind):
#    if CAind>0:
#        LAIind=L*par.sla/CAind
#    else:
#        LAIind=0
#    return LAIind

def lai_tree(L, CAind, sla):
    if CAind>0:
        LAIind=L*sla/CAind
    else:
        LAIind=0
    return LAIind

## missing here:
# actual lai: called in interception() and returns the actual lai of a tree (multiplied with phen)


## old (single PFT)

#def allometry_tree(L,S,H):
#    if ( S<=0.0 or L<=0.0):
#        height=0
#    else:
#        height=S*par.k_latosa/(par.wooddens*L*par.sla)
#
#    if (height>par.height_max):
#        height=par.height_max
#        sm_ind_temp=S
#        S=L*par.height_max*par.wooddens*par.sla/par.k_latosa
#        H=H+sm_ind_temp-S
#
#
#    allometry=par.allom1*np.power(height/par.allom2,par.reinickerp/par.allom3)
#
#    #if allometry > par.CAmax:
#    #    print("CAmax reached")
#
#    CAind=min(allometry,par.CAmax)
#    return S, H, CAind, height



def allometry_tree(L,S,H,sla,CAmax):

    if ( S<=0.0 or L<=0.0):
        height=0
    else:
        height=S*par.k_latosa/(par.wooddens*L*sla)

    if (height>par.height_max):
        height=par.height_max
        sm_ind_temp=S
        S=L*par.height_max*par.wooddens*sla/par.k_latosa
        H=H+sm_ind_temp-S

    allometry=par.allom1*np.power(height/par.allom2,par.reinickerp/par.allom3)

    #if allometry > CAmax:
    #    print("CAmax reached")

    CAind=min(allometry,CAmax)
    return S, H, CAind, height


def fcn(leaf_inc,k1,lm,k3,b,L,H):
    return k1*(b-leaf_inc*lm+H)-np.power((b-leaf_inc*lm)/(L+leaf_inc)*k3,1.0+2.0/par.allom3)
    

def bisect(func, xlow, xhigh, k1,lm,k3,b,L,H, xacc, yacc, maxit):
    ymin=1e09     ## see numeric/bisect.c    
    ylow=func(xlow, k1,lm,k3,b,L,H)
    i=0
    while i < maxit:
        xmid=(xlow+xhigh)*0.5
        if xhigh-xlow<xacc:
            return xmid
        
        ymid=func(xmid,k1,lm,k3,b,L,H)
        if(abs(ymid)<ymin):
            ymin=abs(ymid)
            xmin=xmid
        
        if(abs(ymid)<yacc):
            return xmid
        
        if(ylow*ymid<=0):
            xhigh=xmid
        else:
            xlow=xmid
            ylow=ymid
        i=i+1
    
    return xmin

def leftmostzero(func,x1, x2, k1,lm,k3,b,L,H, xacc, yacc, maxiter):
    if(x2<x1):
        swap=x1
        x1=x2
        x2=swap

    dx=(x2-x1)/par.NSEG
    if(func(x1,k1,lm,k3,b,L,H)<0):
         xmid=x1+dx
         while func(xmid,k1,lm,k3,b,L,H)<0 and xmid<=x2-dx:
              xmid=xmid+dx
    else:  #func>0
        xmid=x1+dx
        while func(xmid,k1,lm,k3,b,L,H)>0 and xmid<=x2-dx:
            xmid=xmid+dx
            
    return bisect(func,xmid-dx,xmid,k1,lm,k3,b,L,H,xacc,yacc,maxiter)



#def allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal):
#    
#    drought1, drought2 = 0, 0
#    lmtorm=par.lmro_ratio*(par.lmro_offset+(1-par.lmro_offset)*min(par.vscal,wscal))
#
#    #bm_inc_ind=bm_inc/N      #### This was feeding in bm_inc from LPJ and divide by N using the state of N in this reduced model.
#        # not used here anymore. It was useful to make the model stable when feeding in the four 
#        # allocated fluxes (model10) because otherwise becomes numerically unstable.
#                              
#    ## alternative: feed in bm_inc_ind directly, i.e. the NPP per individual,
#    ## then the model can decide how many trees there are, and 
#    ## hence also how large the overall NPP per m2 is.
#
#    tinc_H=0
#    tinc_D=0
#    if lmtorm<1.0e-10:
#        print("missing if case")
#    # else
#
#    if height>0:
#        tinc_ind_min_leaf=par.k_latosa*S/(par.wooddens*height*par.sla)-L
#        tinc_ind_min_root=par.k_latosa*S/(par.wooddens*height*par.sla*lmtorm)-R
#    else:
#        tinc_ind_min_leaf=0
#        tinc_ind_min_root=0
#        
#    cmass_deficit=tinc_ind_min_leaf+tinc_ind_min_root-bm_inc_ind
#
#
#    if cmass_deficit>0:
#        cmass_loan=max(min(cmass_deficit*par.CDEBT_MAXLOAN_DEFICIT,S-D)*par.CDEBT_MAXLOAN_MASS,0)
#        bm_inc_ind=bm_inc_ind+cmass_loan
#        tinc_D=cmass_loan
#        
#    else:    # true, except at dry cells!
#        tinc_D=0 
#    
#    
#    if (tinc_ind_min_root>=0 and tinc_ind_min_leaf>=0 and (tinc_ind_min_root+tinc_ind_min_leaf<=bm_inc_ind or bm_inc_ind<=0)):
#        
#        b=S+bm_inc_ind-L/lmtorm+R
#        lm=1+1/lmtorm
#        k1=np.power(par.allom2,2.0/par.allom3)*4.0*(1.0/np.pi)/par.wooddens
#        k3=par.k_latosa/par.wooddens/par.sla
#      
#        x2=(bm_inc_ind-(L/lmtorm-R))/lm
#        if L<1e-10:
#            x1=x2/par.NSEG
#        else:
#            x1=0
#
#        ## bisection:
#        if((x1==0 and x2==0) or b-x1*lm<0 or L+x1<=0 or b-x2*lm<0 or L+x2<=0):
#            tinc_L=0
#        else:
#            tinc_L=leftmostzero(fcn,x1,x2,k1,lm,k3,b,L,H,par.EPSILON_BISECT_x,par.EPSILON_BISECT_y,par.MAXITER_BISECT) # params like in LPJ
#
#        if (tinc_L<0):
#            tinc_R=0
#            
#        else:
#            tinc_R=(tinc_L+L)/lmtorm-R 
#
#        if(bm_inc_ind>0 and tinc_R+tinc_L>bm_inc_ind):
#            tinc_R=bm_inc_ind*tinc_R/(tinc_R+tinc_L)
#            tinc_L=bm_inc_ind*tinc_L/(tinc_R+tinc_L)
#        
#        tinc_S=bm_inc_ind-tinc_L-tinc_R
#        
#
#    else: ##   if NOT: (tinc_ind_min_root>=0 and tinc_ind_min_leaf>=0 and (tinc_ind_min_root+tinc_ind_min_leaf<=bm_inc_ind or bm_inc_ind<=0)):
#
#        ## Abnormal allocation:
#        tinc_L=(bm_inc_ind+R-L/lmtorm)/(1+1/lmtorm)
#        if (tinc_L>0): 
#            ## Sitch 2003:
#            # "In years of stress, the biomass increment may not allow sufficient
#            # allocation to the leaves to fully utilize the current
#            # sapwood (given the constraint implied by Eqn. 1: LA = k_lasa * SA). This
#            # year's production is then allocated to leaves and roots
#            # only, and the excess sapwood mass transferred to the
#            # nonliving heartwood pool."
#            # follows the "pipe model"
#            # "This relationship is based on numerous studies indicating that each unit of leaf area must
#            # be supported by a corresponding area of transport tissue (Shinozaki et al., 1964a, b; Kaufmann & Troendle, 1981;
#            # Waring et al., 1982; Ryan, 1989; Robichaud & Methven, 1992; Berninger & Nikinmaa, 1994)."
#            drought1=1
#            tinc_R=bm_inc_ind-tinc_L
#        else:          
#            ## Sitch 2003:
#            #  "In a year with severe drought
#            # there may be insufficient biomass increment to maintain
#            # both current sapwood mass and leaf mass. In this case all
#            # of the biomass increment is allocated to fine roots and
#            # excess sapwood (...) transferred to the heart-
#            # wood (...) pool."
#            drought2=1
#            tinc_R=bm_inc_ind
#            tinc_L=(R+tinc_R)*lmtorm-L
#            # tinc_L => litter
#        # "excess sapwood" because it is too much for the few leaves,
#        # based on empirical relationship between leaf to sapwood cross-sect area
#        tinc_S=(tinc_L+L)*par.wooddens*height*par.sla/par.k_latosa-S   
#        tinc_H=-tinc_S   # "transferred to the heartwood pool"
#
#    #if(bm_inc_ind>0 and tinc_L>0 and tinc_S>0 and R>0):
#        #print("Computes falloc")
#      
#    #S, H, CAind, height=allometry_tree(L, S, H)
#    #*fpc_inc=fpc_tree(pft)
#    #return isneg_tree(pft)
#
#    #any(x is None for x in tinc_L, tinc_R, tinc_S, tinc_H, tinc_D, drought1, drought2)
#    
##    list1=[tinc_L, tinc_R, tinc_S, tinc_H, tinc_D, drought1, drought2]
##    if None in list1:
##        print("none in allo")
#
#    return tinc_L, tinc_R, tinc_S, tinc_H, tinc_D, drought1, drought2




def allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal, sla):
    
    drought1, drought2 = 0, 0
    lmtorm=par.lmro_ratio*(par.lmro_offset+(1-par.lmro_offset)*min(par.vscal,wscal))

#    #bm_inc_ind=bm_inc/N      #### This was feeding in bm_inc from LPJ and divide by N using the state of N in this reduced model.
#        # not used here anymore. It was useful to make the model stable when feeding in the four 
#        # allocated fluxes (model10) because otherwise becomes numerically unstable.
#                              
#    ## alternative: feed in bm_inc_ind directly, i.e. the NPP per individual,
#    ## then the model can decide how many trees there are, and 
#    ## hence also how large the overall NPP per m2 is.
    
    tinc_H=0
    tinc_D=0
    if lmtorm<1.0e-10:
        print("missing if case")
    # else

    if height>0:
        tinc_ind_min_leaf=par.k_latosa*S/(par.wooddens*height*sla)-L
        tinc_ind_min_root=par.k_latosa*S/(par.wooddens*height*sla*lmtorm)-R
    else:
        tinc_ind_min_leaf=0
        tinc_ind_min_root=0
        
    cmass_deficit=tinc_ind_min_leaf+tinc_ind_min_root-bm_inc_ind


    if cmass_deficit>0:
        cmass_loan=max(min(cmass_deficit*par.CDEBT_MAXLOAN_DEFICIT,S-D)*par.CDEBT_MAXLOAN_MASS,0)
        bm_inc_ind=bm_inc_ind+cmass_loan
        tinc_D=cmass_loan
        
    else:    # true, except at dry cells!
        tinc_D=0 
    
    
    if (tinc_ind_min_root>=0 and tinc_ind_min_leaf>=0 and (tinc_ind_min_root+tinc_ind_min_leaf<=bm_inc_ind or bm_inc_ind<=0)):
        
        b=S+bm_inc_ind-L/lmtorm+R
        lm=1+1/lmtorm
        k1=np.power(par.allom2,2.0/par.allom3)*4.0*(1.0/np.pi)/par.wooddens
        k3=par.k_latosa/par.wooddens/sla
      
        x2=(bm_inc_ind-(L/lmtorm-R))/lm
        if L<1e-10:
            x1=x2/par.NSEG
        else:
            x1=0

        ## bisection:
        if((x1==0 and x2==0) or b-x1*lm<0 or L+x1<=0 or b-x2*lm<0 or L+x2<=0):
            tinc_L=0
        else:
            tinc_L=leftmostzero(fcn,x1,x2,k1,lm,k3,b,L,H,par.EPSILON_BISECT_x,par.EPSILON_BISECT_y,par.MAXITER_BISECT) # params like in LPJ

        if (tinc_L<0):
            tinc_R=0
            
        else:
            tinc_R=(tinc_L+L)/lmtorm-R 

        if(bm_inc_ind>0 and tinc_R+tinc_L>bm_inc_ind):
            tinc_R=bm_inc_ind*tinc_R/(tinc_R+tinc_L)
            tinc_L=bm_inc_ind*tinc_L/(tinc_R+tinc_L)
        
        tinc_S=bm_inc_ind-tinc_L-tinc_R
        

    else: ##   if NOT: (tinc_ind_min_root>=0 and tinc_ind_min_leaf>=0 and (tinc_ind_min_root+tinc_ind_min_leaf<=bm_inc_ind or bm_inc_ind<=0)):

        ## Abnormal allocation:
        tinc_L=(bm_inc_ind+R-L/lmtorm)/(1+1/lmtorm)
        if (tinc_L>0): 
            ## Sitch 2003:
            # "In years of stress, the biomass increment may not allow sufficient
            # allocation to the leaves to fully utilize the current
            # sapwood (given the constraint implied by Eqn. 1: LA = k_lasa * SA). This
            # year's production is then allocated to leaves and roots
            # only, and the excess sapwood mass transferred to the
            # nonliving heartwood pool."
            # follows the "pipe model"
            # "This relationship is based on numerous studies indicating that each unit of leaf area must
            # be supported by a corresponding area of transport tissue (Shinozaki et al., 1964a, b; Kaufmann & Troendle, 1981;
            # Waring et al., 1982; Ryan, 1989; Robichaud & Methven, 1992; Berninger & Nikinmaa, 1994)."
            drought1=1
            tinc_R=bm_inc_ind-tinc_L
        else:          
            ## Sitch 2003:
            #  "In a year with severe drought
            # there may be insufficient biomass increment to maintain
            # both current sapwood mass and leaf mass. In this case all
            # of the biomass increment is allocated to fine roots and
            # excess sapwood (...) transferred to the heart-
            # wood (...) pool."
            drought2=1
            tinc_R=bm_inc_ind
            tinc_L=(R+tinc_R)*lmtorm-L
            # tinc_L => litter
        # "excess sapwood" because it is too much for the few leaves,
        # based on empirical relationship between leaf to sapwood cross-sect area
        tinc_S=(tinc_L+L)*par.wooddens*height*sla/par.k_latosa-S   
        tinc_H=-tinc_S   # "transferred to the heartwood pool"

    #if(bm_inc_ind>0 and tinc_L>0 and tinc_S>0 and R>0):
        #print("Computes falloc")
      
    #S, H, CAind, height=allometry_tree(L, S, H)
    #*fpc_inc=fpc_tree(pft)
    #return isneg_tree(pft)

    #any(x is None for x in tinc_L, tinc_R, tinc_S, tinc_H, tinc_D, drought1, drought2)
    
#    list1=[tinc_L, tinc_R, tinc_S, tinc_H, tinc_D, drought1, drought2]
#    if None in list1:
#        print("none in allo")

    return tinc_L, tinc_R, tinc_S, tinc_H, tinc_D, drought1, drought2





def turnover_tree(L, R, S, H, fL, fR, fS):
    
    turnover_L=-fL*L
    turnover_R=-fR*R
    turnover_S=-fS*S
    turnover_H=fS*S
    turnover=turnover_L+turnover_R+turnover_S
    L=L+turnover_L
    R=R+turnover_R
    S=S+turnover_S
    H=H+turnover_H
    L, R, S, H = crop_pools(L, R, S, H)

    return L, R, S, H, turnover


#def turnover_tree_multiPFT(L, R, S, H, fL, fR, fS):
#
#    NPFT=np.length(L)
#    turnover=np.zeros(NPFT)
#    for PFTind in range(0,NPFT):
#        turnover_L=-fL[PFTind]*L[PFTind]
#        turnover_R=-fR[PFTind]*R[PFTind]
#        turnover_S=-fS[PFTind]*S[PFTind]
#        turnover_H=fS[PFTind]*S[PFTind]
#        turnover[PFTind]=turnover_L+turnover_R+turnover_S
#        L[PFTind]=L[PFTind]+turnover_L
#        R[PFTind]=R[PFTind]+turnover_R
#        S[PFTind]=S[PFTind]+turnover_S
#        H[PFTind]=H[PFTind]+turnover_H
#        L, R, S, H = crop_pools_multiPFT(L, R, S, H)
#
#    return L, R, S, H, turnover


#def tree_mortality(Aperind,N,turnover, L, Nvarfrac):
#    bm_delta=Aperind+turnover
#    if(bm_delta<0):
#        bm_delta=0
#    mort=par.mort_max/(1+par.k_mort*bm_delta/(L*par.sla))
# 
#    # the rest only applies for versions where N is interactive (like in LPJ), not just prescribed
#    if N>-9 and Nvarfrac>-9:
#        Nnew=N*(1-mort)
#        Nmort=Nnew-N
#        N=N+Nmort*Nvarfrac
#        if N<0:
#            N=10**-10
#    
#    return N, mort, bm_delta


## old (single PFT)
#def tree_mortality(Aperind, N, turnover, L, Nvarfrac, mort_prescribed):
#
#    if mort_prescribed<0:
#        bm_delta=Aperind+turnover
#        if(bm_delta<0):
#            bm_delta=0
#        mort=par.mort_max/(1+par.k_mort*bm_delta/(L*par.sla))
#    else:
#        mort=mort_prescribed
#        bm_delta=0
#
#    Nmort=-N*mort
#
#    # the rest only applies for versions where N is interactive (like in LPJ), not just prescribed
#    if N>-9 and Nvarfrac>-9:
#        N=N+Nmort*Nvarfrac
#        if N<0:
#            N=10**-10
#        
#    return N, mort, Nmort, bm_delta



def tree_mortality(Aperind, N, turnover, L, Nvarfrac, mort_prescribed, sla):

    #print(Aperind, N, turnover, L, Nvarfrac, mort_prescribed, sla)
#    if L<=0:
    if mort_prescribed<0:
        bm_delta=Aperind+turnover
        if(bm_delta<0):
            bm_delta=0        

    if L<=0:
        mort=0
        bm_delta=0
    else:
        if mort_prescribed<0:
            mort=par.mort_max/(1+par.k_mort*bm_delta/(L*sla))
        else:
            mort=mort_prescribed
            bm_delta=0
    
    Nmort=-N*mort

    # the rest only applies for versions where N is interactive (like in LPJ), not just prescribed
    if N>-9 and Nvarfrac>-9:
        N=N+Nmort*Nvarfrac
        if N<0:
            N=10**-10
        
    #print(N, mort, Nmort, bm_delta)        
        
        
    return N, mort, Nmort, bm_delta



def tree_establishment(fpc_woody, fpc, n_woody, L, R, S, H, D, N, Nvarfrac, k_est, est_prescribed, estfrac, fpcmax, sla, CAmax, lightextcoeff, Lsapl, Rsapl, Ssapl, Hsapl):
    
## orig LPJ
#    if ( fpc_woody >= par.fpc_tree_max or n_woody <= par.epsilon):
#        S, H, CAind, height = allometry_tree(L,S,H)
#        est=0
#    else:
#        est=k_est*(1-np.exp(-5*(1-fpc_woody)))*(1-fpc_woody)/n_woody
#        est=est*Nvarfrac
#
#        Nold=N
#        Nnew=N+est
#        
#        #print(Nold, est, Nnew)
#        L=(L*Nold+par.Lsapl*est)/Nnew   
#        R=(R*Nold+par.Rsapl*est)/Nnew 
#        S=(S*Nold+par.Ssapl*est)/Nnew 
#        H=(H*Nold+par.Hsapl*est)/Nnew
#        D=D*Nold/Nnew
#        
#        S, H, CAind, height = allometry_tree(L,S,H)        
#        N=Nnew
#        
#        fpc=fpc_tree(L,N,CAind)


    if est_prescribed<par.epsilon:
        # interactive
        if ( fpc_woody >= fpcmax or n_woody <= par.epsilon):
            est=0
        else:
            est=k_est*(1-np.exp(-5*(1-fpc_woody)))*(1-fpc_woody)/n_woody
            est=est*Nvarfrac
    else:
        # prescribed
        est=est_prescribed


    Nold=N
    est_N=est
    Nnew_N=N+est_N
    
    est_pools=est_N*estfrac
    Nnew_pools=N+est_pools
    
    #print(Nold, est, Nnew)
    if Nnew_pools > par.epsilon:    # inserted 27.2.24
        L=(L*Nold+Lsapl*est_pools)/Nnew_pools
        R=(R*Nold+Rsapl*est_pools)/Nnew_pools 
        S=(S*Nold+Ssapl*est_pools)/Nnew_pools 
        H=(H*Nold+Hsapl*est_pools)/Nnew_pools
        D=D*Nold/Nnew_pools
    
    
    S, H, CAind, height = allometry_tree(L,S,H, sla, CAmax)
    
    fpc=fpc_tree(L,Nnew_N,CAind, sla, lightextcoeff)
        
    return Nnew_N, L, R, S, H, D, CAind, height, fpc, est


def tree_establishment_N2(L, R, S, H, D, Nold, est, sla, CAmax, lightextcoeff, Lsapl, Rsapl, Ssapl, Hsapl):
    Nnew=Nold+est
    if Nnew>0 and Nold>=0:
        L=(L*Nold+Lsapl*est)/Nnew
        R=(R*Nold+Rsapl*est)/Nnew 
        S=(S*Nold+Ssapl*est)/Nnew 
        H=(H*Nold+Hsapl*est)/Nnew
        D=D*Nold/Nnew

    S, H, CAind, height = allometry_tree(L,S,H, sla, CAmax)
    
    return L, R, S, H, D, CAind, height




def light_tree(N, fpc, fpcmax, fpc_woody, n_woody, fpc_inc_tree, fpc_inc):

# fpc_inc is PFT-specific, like N and fpc
# fpc_woody, n_woody, fpc_inc_tree are sums over all tree PFTs
# fpcmax is a parameter
    
#  Real ntree;        /* no of tree PFTs currently present */   here_ n_woody
#  Real *fpc_total;   /* total grid FPC for PFTs */
#  Real fpc_inc_tree; /* this years total FPC increment for tree PFTs */
#  Real excess;       /* tree FPC or grass cover to be reduced */
    
    ## from lpj/light.c:
    
    if (fpc_woody>fpcmax):
        if(n_woody):
            f=(fpc_woody-fpcmax)/n_woody
        if(fpc_inc_tree>par.epsilon):
            g=(fpc_woody-fpcmax)/fpc_inc_tree
        
        if (fpc_inc_tree>0.0):  
            excess=g*fpc_inc
        else: 
            excess=f
    
    
    ## from light_tree.c:
    if excess<1e-20:
        N_kill=0
    else:
        N_kill=N*excess/fpc
    
    if(N_kill>N):
        N_kill=N
    
    N=N-N_kill
    return N
    


def adjust_tree(fpc_woody, fpc, L, N, CAind, fpcmax, sla, lightextcoeff):
    i = 0
    if (fpc_woody > fpcmax):

        fpc_end=fpc-(fpc_woody-fpcmax)*fpc/fpc_woody

        while fpc_end < fpc and i < par.MAX_ITER:
            frac=fpc_end/fpc
            N=N*frac
            fpc=fpc_tree(L,N,CAind, sla, lightextcoeff)
            i = i + 1
        
    return N, fpc, i



def crop_pools(L, R, S, H):
    if L<0:
        L=10**-10
    if R<0:
        R=10**-10    
    if S<0:
        S=10**-10    
    if H<0:
        H=10**-10

    return L, R, S, H

#
#def crop_pools_multiPFT(L, R, S, H):
##    if L<0:
##        L=10**-10
##    if R<0:
##        R=10**-10
##    if S<0:
##        S=10**-10
##    if H<0:
##        H=10**-10
#
#    L[L<0]=10**-10
#    R[R<0]=10**-10
#    S[S<0]=10**-10
#    H[H<0]=10**-10
#
#    return L, R, S, H







######################## versions of the reduced model (1 PFT only), used in paper

# interactive allo, interactive dynveg
def LPJ_Cbalance_A1N1(L, R, S, H, D, N, height, bm_inc_ind, wscal, fL, fR, fS, \
                   Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch):
   
    
    #print(par.sla, par.lightextcoeff, par.CAmax)
    
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
    #L, R, S, H = crop_pools(L, R, S, H)
        
    ###### allocation
    AL, AR, AS, AH, AD, drought1, drought2 = allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal, par.sla)
    
    Aperind=AL+AR+AS+AH
    
    L, R, S, H, D = L+AL, R+AR, S+AS, H+AH, D+AD
    L, R, S, H = crop_pools(L, R, S, H)

    S, H, CAind, height=allometry_tree(L,S,H, par.sla, par.CAmax)

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)

    
    ##### mortality
    N, mort, Nmort, bm_delta = tree_mortality(Aperind, N, turnover, L, Nvarfrac, mort_prescribed, par.sla)    

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)

        # diagnosis:
    fpc_beforeest=fpc  # to see how N affects fpc, and how fpc affects est
    N_beforeest=N   # to see how N affects fpc
    
    ######## establishment  (OTHER PFTS MATTER)
    fpc_woody=fpc  # 1 PFT only!
  
    N, L, R, S, H, D, CAind, height, fpc, est = tree_establishment(fpc_woody, fpc, \
                                                n_woody, L, R, S, H, D, N, Nvarfrac, k_est, est_prescribed, estfrac, fpcmax, par.sla, par.CAmax, par.lightextcoeff, par.Lsapl, par.Rsapl, par.Ssapl, par.Hsapl)
    estvsN=est/(N-est)
    
    ####### adjustment (OTHER PFTS MATTER)
    if adj_switch==1:
        fpc_woody=fpc  # 1 PFT only!
        N_old=N
        N, fpc, adj_count = adjust_tree(fpc_woody, fpc, L, N, CAind, fpcmax, par.sla, par.lightextcoeff)
        adjN=(N-N_old)
        adjNrel=(N-N_old)/N_old
    else:
        adjN=0
        adjNrel=0
        adj_count=0

    return L, R, S, H, D, N, fpc, height, AL, AR, AS, AH, drought1, drought2, \
           adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest






# interactive allo, but est and N are read in from LPJ
# sink comes from averaging in saplings, est and Nold are fed in
           # P: as original LPJ, but N is prescribed.

def LPJ_Cbalance_A1N2(L, R, S, H, D, Nold, est, height, bm_inc_ind, wscal, fL, fR, fS):
    
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
        
    ###### allocation    
    AL, AR, AS, AH, AD, drought1, drought2 = allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal, par.sla)
    L, R, S, H, D = L+AL, R+AR, S+AS, H+AH, D+AD
    L, R, S, H = crop_pools(L, R, S, H)
    S, H, CAind, height=allometry_tree(L,S,H, par.sla, par.CAmax)
      
    ######## establishment: rescales individual pools
    L, R, S, H, D, CAind, height = tree_establishment_N2(L, R, S, H, D, Nold, est, par.sla, par.CAmax, par.lightextcoeff, par.Lsapl, par.Rsapl, par.Ssapl, par.Hsapl)

    return L, R, S, H, D, height, AL, AR, AS, AH, drought1, drought2






## only N dynamics, Pools are prescribed
def LPJ_Cbalance_A2N1(L, R, S, H, D, AL, AR, AS, AH, N, fL, fR, fS, \
                   Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch):
    
    ####### turnover
    L_out, R_out, S_out, H_out, turnover = turnover_tree(L, R, S, H, fL, fR, fS) # turnover needed for mortality
    
    ####### allocation (fluxes are read in, but only to compute Aperind which is needed for mortality)
    Aperind=AL+AR+AS+AH
    
    S_out, H_out, CAind, height_out=allometry_tree(L,S,H, par.sla, par.CAmax) # L, S, H from during allo

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)  # L from afterallo
    
    ##### mortality
    N, mort, Nmort, bm_delta = tree_mortality(Aperind, N, turnover, L, Nvarfrac, mort_prescribed, par.sla)

    
    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)

    fpc_beforeest=fpc  # to see how N affects fpc, and how fpc affects est
    N_beforeest=N   # to see how N affects fpc

    ######## establishment
    fpc_woody=fpc  # 1 PFT only!
   
    N, L_out, R_out, S_out, H_out, D_out, CAind, height_out, fpc, est = tree_establishment(fpc_woody, fpc, \
                                                n_woody, L, R, S, H, D, N, Nvarfrac, k_est, est_prescribed, estfrac, fpcmax, par.sla, par.CAmax, par.lightextcoeff, par.Lsapl, par.Rsapl, par.Ssapl, par.Hsapl)
    estvsN=est/(N-est)
    
    
    ####### adjustment
    if adj_switch==1:
        fpc_woody=fpc  # 1 PFT only! else: +fpc2
        N_old=N
        N, fpc, adj_count = adjust_tree(fpc_woody, fpc, L, N, CAind, fpcmax, par.sla, par.lightextcoeff)
        adjN=(N-N_old)
        adjNrel=(N-N_old)/N_old
    elif adj_switch==0:  # no adj is used
        adjN=0
        adjNrel=0
        adj_count=0


    return N, fpc, adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest



## prescibed allo, interactive dynveg
### version F: feed in allocated fluxes but still keep veg dynamics
             ## i.e. with interactive N (similar to model 2.10f)
def LPJ_Cbalance_A0N1(L, R, S, H, N, height, ALperm2, ARperm2, ASperm2, AHperm2, fL, fR, fS, \
                   Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch):
    
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
    L, R, S, H = crop_pools(L, R, S, H)

    ####### allocation (fluxes are read in)

    ## Here, need to stabilise the solution 
    # This works by prescribing fluxes per m2 from LPJ,
    # and then using the actual N from the reduced model to get fluxes per individual again:
    AL=ALperm2/N
    AR=ARperm2/N
    AS=ASperm2/N
    AH=AHperm2/N
    
    Aperind=AL+AR+AS+AH
    
    L, R, S, H = L+AL, R+AR, S+AS, H+AH
    L, R, S, H = crop_pools(L, R, S, H)

    S, H, CAind, height=allometry_tree(L,S,H, par.sla, par.CAmax)

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)
        
    ##### mortality    
    N, mort, Nmort, bm_delta = tree_mortality(Aperind, N, turnover, L, Nvarfrac, mort_prescribed, par.sla)    

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)
    
        # diagnosis
    fpc_beforeest=fpc  # to see how N affects fpc, and how fpc affects est
    N_beforeest=N      # to see how N affects fpc    
    
    ######## establishment
    fpc_woody=fpc       # 1 PFT only!
    N, L, R, S, H, Ddiag, CAind, height, fpc, est = tree_establishment(fpc_woody, fpc, \
                                              n_woody, L, R, S, H, 0, N, Nvarfrac, k_est, est_prescribed, estfrac, fpcmax, par.sla, par.CAmax, par.lightextcoeff, par.Lsapl, par.Rsapl, par.Ssapl, par.Hsapl)

    estvsN=est/(N-est)
    

    ####### adjustment
    if adj_switch==1:
        fpc_woody=fpc  # 1 PFT only! else: +fpc2
        N_old=N
        N, fpc, adj_count = adjust_tree(fpc_woody, fpc, L, N, CAind, fpcmax, par.sla, par.lightextcoeff)
        adjN=(N-N_old)
        adjNrel=(N-N_old)/N_old
    else:
        adjN=0
        adjNrel=0
        adj_count=0


    return L, R, S, H, N, fpc, height, AL, AR, AS, AH, \
           adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest




### Like A0N1 but forced with fluxes per ind
def LPJ_Cbalance_A4N1(L, R, S, H, N, height, AL, AR, AS, AH, fL, fR, fS, \
                   Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch):
    
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
    L, R, S, H = crop_pools(L, R, S, H)

    ####### allocation (fluxes are read in)

    ALperm2=AL*N
    ARperm2=AR*N
    ASperm2=AS*N
    AHperm2=AH*N
        
    
    Aperind=AL+AR+AS+AH
    
    L, R, S, H = L+AL, R+AR, S+AS, H+AH
    L, R, S, H = crop_pools(L, R, S, H)

    S, H, CAind, height=allometry_tree(L,S,H, par.sla, par.CAmax)

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)
        
    ##### mortality    
    N, mort, Nmort, bm_delta = tree_mortality(Aperind, N, turnover, L, Nvarfrac, mort_prescribed, par.sla)    

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)
    
        # diagnosis
    fpc_beforeest=fpc  # to see how N affects fpc, and how fpc affects est
    N_beforeest=N      # to see how N affects fpc    
    
    ######## establishment
    fpc_woody=fpc       # 1 PFT only!
    N, L, R, S, H, Ddiag, CAind, height, fpc, est = tree_establishment(fpc_woody, fpc, \
                                              n_woody, L, R, S, H, 0, N, Nvarfrac, k_est, est_prescribed, estfrac, fpcmax, par.sla, par.CAmax, par.lightextcoeff, par.Lsapl, par.Rsapl, par.Ssapl, par.Hsapl)

    estvsN=est/(N-est)
    

    ####### adjustment
    if adj_switch==1:
        fpc_woody=fpc  # 1 PFT only! else: +fpc2
        N_old=N
        N, fpc, adj_count = adjust_tree(fpc_woody, fpc, L, N, CAind, fpcmax, par.sla, par.lightextcoeff)
        adjN=(N-N_old)
        adjNrel=(N-N_old)/N_old
    else:
        adjN=0
        adjNrel=0
        adj_count=0


    return L, R, S, H, N, fpc, height, ALperm2, ARperm2, ASperm2, AHperm2, \
           adj_count, adjN, adjNrel, mort, Nmort, bm_delta, turnover, est, estvsN, fpc_beforeest, N_beforeest




# interactive allocation, no dynveg; mort is the sink term for pools
    
# Pools are still per ind, but here the sink is that mortality acts like a decay for all pools
   # instead of est that merges in saplings.
# mort is still computed interactively (unless prescribed)
# N can still vary but has no effect on pools per ind
  # for Nvarfract=0, this corresponds to former model 2.11N
    
# N: N is prescribed; mortality affects the pools directly, but establishment does not.
## The old "N"model is like N here with Nvarfrac=0 (no variability)
def LPJ_Cbalance_A1N0(L, R, S, H, D, height, bm_inc_ind, wscal, fL, fR, fS, mort_prescribed):
    
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
    
    ###### allocation
    AL, AR, AS, AH, AD, drought1, drought2 = allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal, par.sla)

    Aperind=AL+AR+AS+AH

    L, R, S, H, D = L+AL, R+AR, S+AS, H+AH, D+AD
    L, R, S, H = crop_pools(L, R, S, H)
    S, H, CAind, height=allometry_tree(L,S,H, par.sla, par.CAmax)
    
    if mort_prescribed<0:
        Ndiag, mort, Nmortdiag, bm_delta = tree_mortality(Aperind, -999, turnover, L, -999, mort_prescribed, par.sla)
    else:
        mort=mort_prescribed
        bm_delta=0
    

    L=L-mort*L
    R=R-mort*R
    S=S-mort*S
    H=H-mort*H
    
    return L, R, S, H, D, height, AL, AR, AS, AH, drought1, drought2, mort, bm_delta, turnover


##### input: A perm2
## no allo, no dynveg
## This is like model 2.10N, but again, N can vary, but this just scales everything
#def LPJ_Cbalance_A0N0(L, R, S, H, N, ALperm2, ARperm2, ASperm2, AHperm2, fL, fR, fS, mort_prescribed):
#        
#    ####### turnover
#    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
#    
#    AL=ALperm2/N
#    AR=ALperm2/N
#    AS=ALperm2/N
#    AH=ALperm2/N
#    
#    ###### allocation
#    L, R, S, H = L+AL, R+AR, S+AS, H+AH
#    L, R, S, H = crop_pools(L, R, S, H)
#    S, H, CAind, height=allometry_tree(L,S,H)
#    
#    ##### mortality
#    #Aperm2=(AL+AR+AS+AH)*N
#    Aperind=AL+AR+AS+AH
#
#    if mort_prescribed<0:
#        Ndiag, mort, Nmortdiag, bm_delta = tree_mortality(Aperind, -999, turnover, L, -999, mort_prescribed)    
#    else:
#        mort=mort_prescribed
#        bm_delta=0
#        
#    L=L-mort*L
#    R=R-mort*R
#    S=S-mort*S
#    H=H-mort*H
#    
#    return L, R, S, H, mort, bm_delta, turnover
#


### input: A per ind
# no allo, no dynveg
# This is like model 2.10N, but again, N can vary, which scales the pools from per ind to per m2
def LPJ_Cbalance_A0N0(L, R, S, H, AL, AR, AS, AH, fL, fR, fS, mort_prescribed):
        
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
    
    ###### allocation
    L, R, S, H = L+AL, R+AR, S+AS, H+AH
    L, R, S, H = crop_pools(L, R, S, H)
    S, H, CAind, height=allometry_tree(L,S,H, par.sla, par.CAmax)
    
    ##### mortality
    #Aperm2=(AL+AR+AS+AH)*N
    Aperind=AL+AR+AS+AH

    if mort_prescribed<0:
        Ndiag, mort, Nmortdiag, bm_delta = tree_mortality(Aperind, -999, turnover, L, -999, mort_prescribed, par.sla)
    else:
        mort=mort_prescribed
        bm_delta=0
        
    L=L-mort*L       
    R=R-mort*R
    S=S-mort*S
    H=H-mort*H
    
    return L, R, S, H, mort, bm_delta, turnover






################### reduced versions for multiple PFTs

def LPJ_Cbalance_A1N1_multiPFT_1dimtest(L, R, S, H, D, N, height, bm_inc_ind, wscal, fpc_woody_beforeest, fpc_rest_beforeest, fpc_woody_afterest, fpc_rest_afterest, fL, fR, fS, \
                   Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch):
   
    ## for each tree PFT...
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H, fL, fR, fS)
    
    ###### allocation
    AL, AR, AS, AH, AD, drought1, drought2 = allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal, par.sla)
    
    Aperind=AL+AR+AS+AH
    
    L, R, S, H, D = L+AL, R+AR, S+AS, H+AH, D+AD
    L, R, S, H = crop_pools(L, R, S, H)

    S, H, CAind, height=allometry_tree(L,S,H, par.sla, par.CAmax)

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)

    
    ##### mortality
    N, mort, Nmort, bm_delta = tree_mortality(Aperind, N, turnover, L, Nvarfrac, mort_prescribed, par.sla)    

    fpc=fpc_tree(L,N,CAind, par.sla, par.lightextcoeff)

    #fpc_beforeest=fpc  # to see how N affects fpc, and how fpc affects est
    #N_beforeest=N   # to see how N affects fpc
    
    ### end
    
    
    ######## establishment  (OTHER PFTS MATTER)
    #fpc_woody=fpc  # 1 PFT only!
  
    N, L, R, S, H, D, CAind, height, fpc, est = tree_establishment(fpc_woody_beforeest, fpc, \
                                                n_woody, L, R, S, H, D, N, Nvarfrac, k_est, est_prescribed, estfrac, fpcmax, par.sla, par.CAmax, par.lightextcoeff, par.Lsapl, par.Rsapl, par.Ssapl, par.Hsapl)

    estvsN=est/(N-est)
    
    ####### adjustment (OTHER PFTS MATTER)
    if adj_switch==1:
        #fpc_woody=fpc  # 1 PFT only!
        N_old=N
        N, fpc, adj_count = adjust_tree(fpc_woody_afterest, fpc, L, N, CAind, fpcmax, par.sla, par.lightextcoeff)
        adjN=(N-N_old)
        adjNrel=(N-N_old)/N_old
    else:
        adjN=0
        adjNrel=0
        adj_count=0


    return L, R, S, H, D, N, fpc, height, adj_count, adjN, adjNrel, mort, Nmort, est, estvsN




def LPJ_Cbalance_A1N1_multiPFT(L, R, S, H, D, N, height, bm_inc_ind, wscal, fL, fR, fS, \
                               Nvarfrac, k_est, est_prescribed, estfrac, n_woody, mort_prescribed, fpcmax, adj_switch, NPFT):
       
    ###L, R, S, H, D, N, height, bm_inc_ind, wscal, fL, fR, fS, are vectors with length NPFT
    ## the rest are scalar parameters
    fpc=np.zeros(NPFT)
    
    #turnover=np.zeros(NPFT)
    mort=np.zeros(NPFT)
    Nmort=np.zeros(NPFT)
    
    #print('before allo: ', L)

    
    for PFTind in range(0,NPFT):  # NPFT is the number of PFTs, not the index of the last one.

        if bm_inc_ind[PFTind]>0 or wscal[PFTind]>0: # check if PFTind is active
    
            ####### turnover
            L[PFTind], R[PFTind], S[PFTind], H[PFTind], turnover = turnover_tree(L[PFTind], R[PFTind], S[PFTind], H[PFTind], fL[PFTind], fR[PFTind], fS[PFTind])
            #L, R, S, H, turnover = turnover_tree_multiPFT(L, R, S, H, fL, fR, fS, NPFT)
    
            ###### allocation
            AL, AR, AS, AH, AD, drought1, drought2 = allocation_tree(L[PFTind], R[PFTind], S[PFTind], H[PFTind], D[PFTind], height[PFTind], bm_inc_ind[PFTind], wscal[PFTind], par.sla[PFTind])
            
            Aperind=AL+AR+AS+AH
            
            L[PFTind], R[PFTind], S[PFTind], H[PFTind], D[PFTind] = L[PFTind]+AL, R[PFTind]+AR, S[PFTind]+AS, H[PFTind]+AH, D[PFTind]+AD
    
            L[PFTind], R[PFTind], S[PFTind], H[PFTind] = crop_pools(L[PFTind], R[PFTind], S[PFTind], H[PFTind])

            S[PFTind], H[PFTind], CAind, height[PFTind]=allometry_tree(L[PFTind],S[PFTind],H[PFTind], par.sla[PFTind], par.CAmax[PFTind])
        
            fpc[PFTind]=fpc_tree(L[PFTind],N[PFTind],CAind, par.sla[PFTind], par.lightextcoeff[PFTind])
                    
            ##### mortality
            N[PFTind], mort[PFTind], Nmort[PFTind], bm_delta = tree_mortality(Aperind, N[PFTind], turnover, L[PFTind], Nvarfrac, mort_prescribed, par.sla[PFTind])
        
            fpc[PFTind]=fpc_tree(L[PFTind],N[PFTind],CAind, par.sla[PFTind], par.lightextcoeff[PFTind])
        
            #fpc_beforeest=fpc  # to see how N affects fpc, and how fpc affects est
            #N_beforeest=N   # to see how N affects fpc
            
        else:
            L[PFTind], R[PFTind], S[PFTind], H[PFTind], D[PFTind], \
            N[PFTind], mort[PFTind], Nmort[PFTind], fpc[PFTind] = np.zeros(9)

    ### end
    
    #print('before est: ', L)

    
    ######## establishment  (OTHER PFTS MATTER)

    fpc_woody=np.sum(fpc)  # 1 PFT only!
    est=np.zeros(NPFT)
    estvsN=np.zeros(NPFT)

    for PFTind in range(0,NPFT):

        #print('input: ', PFTind, bm_inc_ind[PFTind], wscal[PFTind])
        
        #if bm_inc_ind[PFTind]>0 or wscal[PFTind]>0: # check if PFTind is active
            
        N[PFTind], L[PFTind], R[PFTind], S[PFTind], H[PFTind], D[PFTind], CAind, \
        height[PFTind], fpc[PFTind], est[PFTind] \
            = tree_establishment(fpc_woody, fpc[PFTind], \
                                 n_woody, L[PFTind], R[PFTind], S[PFTind], H[PFTind], D[PFTind], \
                                 N[PFTind], Nvarfrac, k_est, est_prescribed, estfrac, fpcmax, \
                                 par.sla[PFTind], par.CAmax[PFTind], par.lightextcoeff[PFTind], \
                                 par.Lsapl[PFTind], par.Rsapl[PFTind], par.Ssapl[PFTind], par.Hsapl[PFTind])
        estvsN[PFTind]=est[PFTind]/(N[PFTind]-est[PFTind])
    
        #print('L after est: ', L[PFTind])
        
#        else:
#            L[PFTind], R[PFTind], S[PFTind], H[PFTind], D[PFTind], \
#            N[PFTind], mort[PFTind], Nmort[PFTind], fpc[PFTind] = np.zeros(9)        




    #print('before adj: ', L)
    
    ####### adjustment (OTHER PFTS MATTER)
    adj_count=np.zeros(NPFT)
    adjN=np.zeros(NPFT)
    adjNrel=np.zeros(NPFT)

    if adj_switch==1:
        fpc_woody=np.sum(fpc) # scalar
        #fpc_woody=fpc  # 1 PFT only!
        N_old=N # vector
        
        for PFTind in range(0,NPFT):
            #if bm_inc_ind[PFTind]>0 or wscal[PFTind]>0: # check if PFTind is active
            N[PFTind], fpc[PFTind], adj_count[PFTind] = adjust_tree(fpc_woody, fpc[PFTind], \
             L[PFTind], N[PFTind], CAind, fpcmax, par.sla[PFTind], par.lightextcoeff[PFTind])
            adjN[PFTind]=(N[PFTind]-N_old[PFTind])

            adjNrel[PFTind]=(N[PFTind]-N_old[PFTind])/N_old[PFTind]

        
#            else:
#                L[PFTind], R[PFTind], S[PFTind], H[PFTind], D[PFTind], \
#                N[PFTind], mort[PFTind], Nmort[PFTind], fpc[PFTind] = np.zeros(9)
#                adjN[PFTind]=0
#                adjNrel[PFTind]=0
        
    else:
        adjN=0
        adjNrel=0
        adj_count=0

    #print('after adj: ', L)


    ## all output variables are vectors
    return L, R, S, H, D, N, fpc, height, adj_count, adjN, adjNrel, mort, Nmort, est, estvsN

