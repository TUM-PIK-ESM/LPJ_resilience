#!/bin/ksh

## here, will cut off the first year from all fluxes. Otherwise cannot be compared to pool dynamics!

## the years are exp specific!
#explist="sba0034 sba0035 sba0036"
##year_ini=2979     # already assumes that first 1000 years were cut off!
##year_fin=11978

#explist="sba0041"
#year_ini=1979
#year_fin=1985

#explist="sba0042"
#explist="sba0044"
#explist="sba0042 sba0044 sba0047 sba0048 sba0049 sba0052 sba0053 sba0054"

#explist="sba0042 sba0044 sba0047 sba0048 sba0052 sba0053 sba0054"

# sba0042: all PFTs
# sba0044: PFT1
# sba0047: PFT2
# sba0048: PFT3

#explist="sba0042 sba0044 sba0053 sba0056"
#explist="sba0056"
#explist="sba0056 sba0061 sba0062 sba0063 sba0064 sba0065 sba0066 sba0067 sba0068"
explist="sba0082 sba0084 sba0087"

year_ini=1979
year_fin=11978

(( yinip1 = ${year_ini} + 1 ))
(( yfinm1 = ${year_fin} - 1 ))

for exp in ${explist}; do


 if [[ ${exp} == "sba0084" || ${exp} == "sba0087" ]]; then
   PFTmax=1
 else
   PFTmax=8
 fi
 (( PFTmaxmin1 = ${PFTmax} - 1 ))


cd /home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/${exp}

#rm -f n_woody.nc
if [[ ! -f n_woody.nc ]]; then
  cdo vertmax pft_n_est.nc n_woody.nc
fi


#rm -f n_woody_proxy.nc
if [[ ! -f n_woody_proxy.nc ]]; then
# ## count active woody types
  cdo sellevel,1/${PFTmax} pft_nind_timmean.nc pft_nind_timmean_1-${PFTmax}.nc
##  cdo gtc,1e-20 pft_nind_timmean_1-8.nc pft_woody.nc 
  cdo gtc,1e-10 pft_nind_timmean_1-${PFTmax}.nc pft_woody.nc  # original
##  cdo gtc,1e-5 pft_nind_timmean_1-8.nc pft_woody.nc 
##  cdo gtc,1e-3 pft_nind_timmean_1-8.nc pft_woody.nc
  cdo vertsum pft_woody.nc n_woody_proxy.nc
  rm pft_nind_timmean_1-${PFTmax}.nc
fi



#rm -f fpc_woody.nc
if [[ ! -f fpc_woody_timmean.nc ]]; then
  cdo sellevel,1/${PFTmax} fpc.nc fpc_1-${PFTmax}.nc
  cdo mul fpc_1-${PFTmax}.nc pft_woody.nc pft_fpc_woody.nc
  cdo vertsum pft_fpc_woody.nc fpc_woody.nc
  rm pft_fpc_woody.nc fpc_1-${PFTmax}.nc
  cdo timmean fpc_woody.nc fpc_woody_timmean.nc
fi

#### total individuals per m2
if [[ ! -f nind.nc ]]; then
  cdo vertsum pft_nind.nc nind.nc
  cdo timmean nind.nc nind_timmean.nc
  cdo timmean pft_nind.nc pft_nind_timmean.nc
fi


### carbon pools per m2: g/ind x ind/m2
for pool in cleaf csapw croot chawo; do
 if [[ ! -f ${pool}perm2.nc ]]; then
   cdo mul pft_nind.nc pft_${pool}.nc pft_${pool}perm2.nc
   cdo timmean pft_${pool}perm2.nc pft_${pool}perm2_timmean.nc
   cdo vertsum pft_${pool}perm2.nc ${pool}perm2.nc
 fi
done

# total C perm2 = vegc!

## veg carbon per ind, PFT-specific
if [[ ! -f pft_vegcperind_timmean.nc ]]; then
  cdo div pft_vegc.nc pft_nind.nc pft_vegcperind.nc
  cdo timmean pft_vegcperind.nc pft_vegcperind_timmean.nc
fi

## veg carbon per ind, total
if [[ ! -f vegcperind_timmean.nc ]]; then
  cdo div vegc.nc nind.nc vegcperind.nc
  cdo timmean vegcperind.nc vegcperind_timmean.nc
fi

# veg carbon per m2, total: vegc
# veg carbon per m2, per PFT: pft_vegc

### fraction of vegetation carbon in each pft
#if [[ ! -f pft_vegcfrac.nc ]]; then
#  cdo div pft_vegc.nc vegc.nc pft_vegcfrac.nc
#fi






#
### fractionation of pools
if [[ ! -f pft_SH_frac_timmean.nc ]]; then

 for pool in cleaf croot csapw chawo; do
  cdo div pft_${pool}_timmean.nc pft_vegcperind_timmean.nc pft_${pool}_frac_timmean.nc
 done

  cdo add pft_csapw_timmean.nc pft_chawo_timmean.nc pft_SH.nc
#  ##cdo add pft_cleaf.nc pft_croot.nc pft_LR.nc
#  ##cdo add pft_LR.nc pft_SH.nc pft_C.nc

  cdo div pft_SH.nc pft_vegcperind_timmean.nc pft_SH_frac_timmean.nc
  rm pft_SH.nc
  
fi


# overall fractionation (not PFT-specific):





#
## total est
#if [[ ! -f est.nc ]]; then
#  cdo vertsum pft_est.nc est.nc
#fi
#
## net mortality over all PFTs;
#if [[ ! -f mort.nc ]]; then
#  cdo vertsum pft_nind_beforemort.nc nind_beforemort.nc 
#  cdo vertsum pft_nind_aftermort.nc nind_aftermort.nc  
#  cdo div nind_aftermort.nc nind_beforemort.nc ratio.nc
#  cdo addc,1 -mulc,-1 ratio.nc mort.nc
#  rm ratio.nc
#  cdo setmisstoc,0 mort.nc mort0.nc
#  mv mort0.nc mort.nc
#fi
#
#
#
######## adj diagnosis
#if [[ ! -f pft_adjNrel_timmean.nc ]]; then
#
#  # wird gebraucht: nind before and after adj
#  #    adjN=-(N-N_old)
#  #    adjNrel=-(N-N_old)/N_old
#
#  cdo sub pft_nind_beforeadj.nc pft_nind_afteradj.nc pft_adjN.nc
#
#  cdo div pft_adjN.nc pft_nind_beforeadj.nc pft_adjNrel.nc
#
#  cdo setname,pft_nind pft_adjNrel.nc pft_adjNrel_named.nc
#  mv pft_adjNrel_named.nc pft_adjNrel.nc
#  cdo setname,pft_nind pft_adjN.nc pft_adjN_named.nc
#  mv pft_adjN_named.nc pft_adjN.nc
#
#  cdo timmean pft_adjN.nc pft_adjN_timmean.nc
#  cdo timmean pft_adjNrel.nc pft_adjNrel_timmean.nc
#
#fi
#
## now again for total N:
#if [[ ! -f adjN_timmean.nc ]]; then
#  cdo vertsum pft_adjN.nc adjN.nc
#  
## cdo setctomiss,0 pft_nind_aftermort.nc miss1.nc
#  
#  cdo vertsum pft_adjN_timmean.nc adjN_timmean.nc
#fi
#if [[ ! -f adjNrel_timmean.nc ]]; then
#  cdo vertsum pft_nind_beforeadj.nc nind_beforeadj.nc
#  cdo div adjN.nc nind_beforeadj.nc adjNrel.nc
#  cdo timmean adjNrel.nc adjNrel_timmean.nc
#fi
#
#
### correlations
#if [[ ! -f pft_Cest_timcorr.nc ]]; then
#  cdo timcor pft_nind.nc pft_vegcperind.nc pft_CN_timcorr.nc
#  cdo timcor pft_nind.nc pft_vegc.nc pft_Cperm2N_timcorr.nc
#  #cdo sellevel,1 fpc_afterest.nc fpc_afterest_lev1.nc
#  #cdo sellevel,1 pft_nind_afterest.nc pft_nind_afterest_lev1.nc
#  #cdo timcor pft_nind_afterest_lev1.nc fpc_afterest_lev1.nc pft_Nfpc_timcorr.nc
#  cdo timcor pft_nind.nc fpc.nc pft_Nfpc_timcorr.nc
#  cdo timcor pft_nind_beforeest.nc pft_est.nc pft_Nest_timcorr.nc 
#  cdo timcor pft_est.nc pft_vegcperind.nc pft_Cest_timcorr.nc
#  rm fpc_afterest_lev1.nc pft_nind_afterest_lev1.nc
#fi
#
#
#
#### total fpc
##if [[ ! -f fpc_timmean.nc ]]; then
##  cdo vertsum fpc.nc fpc_vertsum.nc
##  cdo timmean fpc_vertsum.nc fpc_vertsum_timmean.nc
##  cdo timmean fpc.nc fpc_timmean.nc
##fi
#
#
### diagnose bm_inc_ind
##if [[ ! -f pft_bm_inc_ind.nc ]]; then
##  cdo div pft_bm_inc_carbon_beforeallo.nc pft_nind_beforemort.nc pft_bm_inc_ind.nc
##fi
#
#### diversity in terms of fpc and nind (Shannon entropy)
# ## S = minus sum of pi * ln (pi)
## first need to compute fractions:
##for var in fpc pft_nind; do
##if [[ ! -f Shannon_${var}.nc ]]; then
## cdo div ${var}_timmean.nc ${var}_vertsum_timmean.nc ${var}_rel.nc
## cdo setrtomiss,-1,1e-20 ${var}_rel.nc ${var}_rel_miss.nc
## cdo ln ${var}_rel_miss.nc ${var}_ln.nc
## cdo mul ${var}_rel.nc ${var}_ln.nc mul.nc
## cdo vertsum mul.nc mul2.nc
## cdo mulc,-1 mul2.nc Shannon_${var}.nc
## rm ${var}_rel.nc ${var}_rel_miss.nc ${var}_ln.nc mul.nc mul2.nc
##fi
##done
#
#
#### check.  sum of allocated fluxes is not identical to bm_inc because of some excess, debt, legacy stuff. but very small differences.
###cdo enssum pft_tinc_ind_cleaf.nc pft_tinc_ind_croot.nc pft_tinc_ind_csapw.nc pft_tinc_ind_chawo.nc A_LPJ.nc
##cdo enssum pft_tinc_ind_cleaf.nc pft_tinc_ind_croot.nc pft_tinc_ind_csapw.nc pft_tinc_ind_chawo.nc pft_tinc_ind_cdebt.nc A_LPJ.nc
##cdo sub A_LPJ.nc bm_inc_ind.nc test.nc
##cdo timmean test.nc test_timmean.nc
#


done #exp

exit
