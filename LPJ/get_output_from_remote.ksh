#!/bin/ksh

# fpc hat 12 levels, 7-9 und 12 sind zero!
# pft_npp und pft_vegc haben 43 levels, aber ab 11 ist alles zero.

calc_new=0


#explist="sba0034 sba0035 sba0036"
#year_ini=2979     # already assumes that first 1000 years were cut off!
#year_fin=11978

#explist="sba0041"
#year_ini=1979
#year_fin=1985

## PFT coexistence checks
#explist="sba0042 sba0044 sba0047 sba0048 sba0049 sba0052 sba0053 sba0054"
#explist="sba0053"


#### get all vars I need:
#vars="npp gpp vegc agb litc pft_cleaf pft_csapw pft_croot pft_chawo rh soilc flux_estab fpc pft_mort pft_vegc nv_lai pft_laimax pft_npp pft_nind soilc_slow litfallc litfallc_wood rh_litter soilc_layer cshift_fast_nv cshift_slow_nv tau_cleaf tau_croot tau_csapw pft_turnover_ind pft_bm_inc_carbon pft_bm_delta pft_nind_frac pft_est pft_nind_afterest pft_nind_aftermort pft_nind_afterlight pft_nind_afteradj pft_nind_afterred pft_nind_aftermix pft_nind_beforeest pft_nind_beforelight pft_nind_beforeadj pft_nind_beforered pft_nind_beforemix pft_count_est pft_count_mort pft_count_adj pft_est_cum pft_est_nind pft_fpc_type pft_n_est pft_cleaf_afterest pft_croot_afterest pft_csapw_afterest pft_chawo_afterest pft_cleaf_beforeest pft_croot_beforeest pft_csapw_beforeest pft_chawo_beforeest pft_sla pft_mort_max pft_cleaf_duringmort pft_nind_aftermort pft_count_mort  pft_cleaf_beforeturn pft_cleaf_afterturn fpc_beforemort fpc_aftermort fpc_beforeest fpc_afterest fpc_beforeallo fpc_afterallo fpc_beforeadj fpc_afteradj fpc_beforelight fpc_afterlight pft_CAind_beforeest pft_CAind_afterest"
#vars="pft_tinc_ind_cleaf pft_tinc_ind_croot pft_tinc_ind_csapw pft_tinc_ind_chawo pft_nind_beforemort pft_cleaf_beforeallo pft_cleaf_afterallo pft_CAind_beforeallo pft_CAind_afterallo pft_croot_beforeallo pft_csapw_beforeallo pft_chawo_beforeallo pft_cdebt_beforeallo pft_height pft_wscal excessc_beforeallo excessc_afterallo pft_bm_inc_carbon_beforeallo pft_bm_inc_carbon_afterallo pft_tinc_ind_cdebt"

#vars="pft_wscal"


#vars="npp pft_npp"


#explist="sba0046"
#year_ini=1979
#year_fin=2978
#vars="agb fpc"


#explist="sba0056 sba0061 sba0062 sba0063 sba0064 sba0065 sba0066 sba0067 sba0068 sba0042 sba0044 sba0047 sba0048 sba0049 sba0052 sba0053 sba0054"


## mean state
#vars="fpc pft_wscal pft_nind_beforemort pft_bm_inc_carbon_beforeallo pft_nind pft_bm_inc_carbon n_woody vegc agb     pft_nind_beforelight pft_nind_afterlight pft_nind_afteradj pft_nind_afterest pft_nind_aftermix pft_nind_aftermort pft_nind_afterred pft_nind_beforeadj pft_nind_beforeest pft_nind_beforemix pft_nind_beforemort pft_nind_beforered phen_water pft_n_est pft_chawo pft_csapw pft_cleaf pft_croot"

### shorter period
#explist="sba0042"
#vars="vegc500yrs"


#vars="vegc pft_nind pft_chawo pft_csapw pft_cleaf pft_croot pft_cdebt pft_height"


## mesosouthamerica:
#explist="sba0080"
#vars="vegc agb fpc"

# Nian 24 hypothese
#explist="sba0046"
#vars="agb_tree"


## sba0082-87
#explist="sba0082 sba0084 sba0087"
#explist="sba0080"
#vars="agb_tree agb vegc pft_nind pft_n_est fpc pft_vegc pft_chawo pft_cleaf pft_croot pft_csapw pft_bm_delta pft_wscal pft_bm_inc_carbon_afterallo pft_bm_inc_carbon_beforeallo pft_nind_beforemort pft_nind_beforeest pft_nind_afterest pft_est pft_mort pft_tinc_ind_cleaf pft_tinc_ind_croot pft_tinc_ind_csapw pft_tinc_ind_chawo pft_mort"


## runs with input2
#explist="sba0056 sba0061 sba0062 sba0063 sba0064 sba0065 sba0066 sba0067 sba0068"
#explist="sba0042 sba0044"
#vars="pft_nind_beforemort pft_bm_inc_carbon_beforeallo pft_n_est" 
#vars="agb pft_wscal swc1 swc2 swc pft_vegc pft_chawo pft_cleaf pft_croot pft_csapw pft_npp fpc"
#vars="pft_wscal"


#vars="lai_eff gpp npp fpc_woody" ## to check if resilience also works for other variables



#explist="sba0070 sba0071 sba0072 sba0073 sba0074 sba0075 sba0076 sba0077 sba0078"
#vars="gpp pft_nind agb vegc pft_chawo pft_cleaf pft_croot pft_csapw pft_npp fpc"




## check agb_tree in other runs: => this is not the reason for the shifted MAPcrit
#explist="sba0042"
#vars="agb_tree"




#explist="sba0044 sba0047"


#vars="pft_bm_inc_carbon_beforeallo pft_nind pft_wscal pft_nind_beforemort"


## second: run combine_variables to get pft_bm_inc_ind

## third: run this
#vars="pft_bm_inc_ind pft_nind pft_wscal"


# sba0042: all PFTs
# sba0044: PFT1
# sba0047: PFT2
# sba0048: PFT3

## AC of vegc, and vegc per ind, and nind; each in total and per PFT:
#explist="sba0042 sba0044 sba0047 sba0048 sba0052 sba0053 sba0054"
##vars="pft_vegc pft_nind vegc nind vegcperind pft_vegcperind"
##vars="vegc nind agb"

####

for var in ${vars}; do
  for exp in ${explist}; do

    mkdir -p ${exp}
    cd ${exp}
 
    if [[ ! -f ${var}.nc ]]; then
      scp bathiany@cluster.pik-potsdam.de:/p/projects/tipes/bathiany/LPJmL/${exp}/output/${var}.nc .

      if [[ ( ${exp} == sba0034 || ${exp} == sba0035 || ${exp} == sba0036 || ${exp} == sba0041 ) && -f ${var}.nc ]]; then 
        cdo selyear,${year_ini}/${year_fin} ${var}.nc ${var}_cut.nc
        mv ${var}_cut.nc ${var}.nc
      fi
      
      if [[ ${var} == "fpc" ]]; then
        cdo sellevel,2/12 ${var}.nc ${var}_cut.nc   # there are 11 PFTs, but first level in this file is natural fraction per veg fraction. result will be called level 1-11 after download
        mv ${var}_cut.nc ${var}.nc
      elif [[ ! ${var} == "swc" ]]; then
        cdo sellevel,0/11 ${var}.nc ${var}_cut.nc   # from level 0 because files without levels have only 0
        mv ${var}_cut.nc ${var}.nc
      fi
    fi

    cd ..
  done #var
done #exp

exit
