#!/bin/ksh

# fpc hat 12 levels, 7-9 und 12 sind zero!
# pft_npp und pft_vegc haben 43 levels, aber ab 11 ist alles zero.

calc_new=0

res=y  # always for LPJ because pools and N have yearly steps

#explist="sba0013 sba0015 sba0017 sba0021 sba0016 sba0022 sba0023 sba0033 sba0034"

#explist="sba0034 sba0035 sba0036"
#year_ini=2979     # already assumes that first 1000 years were cut off!
#year_fin=11978

#explist="sba0041"
#year_ini=1979
#year_fin=1985

## PFT coexistence checks
#explist="sba0042 sba0044 sba0047 sba0048 sba0049 sba0052 sba0053 sba0054"
#explist="sba0053"
#year_ini=1979
#year_fin=11978

#### get all vars I need:
#vars="npp gpp vegc agb litc pft_cleaf pft_csapw pft_croot pft_chawo rh soilc flux_estab fpc pft_mort pft_vegc nv_lai pft_laimax pft_npp pft_nind soilc_slow litfallc litfallc_wood rh_litter soilc_layer cshift_fast_nv cshift_slow_nv tau_cleaf tau_croot tau_csapw pft_turnover_ind pft_bm_inc_carbon pft_bm_delta pft_nind_frac pft_est pft_nind_afterest pft_nind_aftermort pft_nind_afterlight pft_nind_afteradj pft_nind_afterred pft_nind_aftermix pft_nind_beforeest pft_nind_beforelight pft_nind_beforeadj pft_nind_beforered pft_nind_beforemix pft_count_est pft_count_mort pft_count_adj pft_est_cum pft_est_nind pft_fpc_type pft_n_est pft_cleaf_afterest pft_croot_afterest pft_csapw_afterest pft_chawo_afterest pft_cleaf_beforeest pft_croot_beforeest pft_csapw_beforeest pft_chawo_beforeest pft_sla pft_mort_max pft_cleaf_duringmort pft_nind_aftermort pft_count_mort  pft_cleaf_beforeturn pft_cleaf_afterturn fpc_beforemort fpc_aftermort fpc_beforeest fpc_afterest fpc_beforeallo fpc_afterallo fpc_beforeadj fpc_afteradj fpc_beforelight fpc_afterlight pft_CAind_beforeest pft_CAind_afterest"
#vars="pft_tinc_ind_cleaf pft_tinc_ind_croot pft_tinc_ind_csapw pft_tinc_ind_chawo pft_nind_beforemort pft_cleaf_beforeallo pft_cleaf_afterallo pft_CAind_beforeallo pft_CAind_afterallo pft_croot_beforeallo pft_csapw_beforeallo pft_chawo_beforeallo pft_cdebt_beforeallo pft_height pft_wscal excessc_beforeallo excessc_afterallo pft_bm_inc_carbon_beforeallo pft_bm_inc_carbon_afterallo pft_tinc_ind_cdebt"

#vars="pft_wscal"
#typelist="timmean"

#vars="npp pft_npp"
#typelist="timstd"

#explist="sba0046"
#year_ini=1979
#year_fin=2978
#vars="agb fpc"


#explist="sba0056 sba0061 sba0062 sba0063 sba0064 sba0065 sba0066 sba0067 sba0068 sba0042 sba0044 sba0047 sba0048 sba0049 sba0052 sba0053 sba0054"
#year_ini=1979
#year_fin=11978

## mean state
#vars="fpc pft_wscal pft_nind_beforemort pft_bm_inc_carbon_beforeallo pft_nind pft_bm_inc_carbon n_woody vegc agb     pft_nind_beforelight pft_nind_afterlight pft_nind_afteradj pft_nind_afterest pft_nind_aftermix pft_nind_aftermort pft_nind_afterred pft_nind_beforeadj pft_nind_beforeest pft_nind_beforemix pft_nind_beforemort pft_nind_beforered phen_water pft_n_est pft_chawo pft_csapw pft_cleaf pft_croot"

#vars="vegc"
##typelist="timmean AC${res}L10"
#typelist="timmin timmean"


### shorter period
#explist="sba0042"
#vars="vegc500yrs"
#typelist="timmean AC${res}L01 AC${res}L10"
#typelist="timmin timmean"


#vars="vegc pft_nind pft_chawo pft_csapw pft_cleaf pft_croot pft_cdebt pft_height"
#typelist="timmean"

#explist="sba0056"
#vars="pft_est"
#typelist="timmean"

## recovery rate
#explist="sba0070"
explist="sba0070 sba0071 sba0072 sba0073 sba0074 sba0075 sba0076 sba0077 sba0078"
#vars="gpp pft_nind agb vegc pft_chawo pft_cleaf pft_croot pft_csapw pft_npp fpc"
vars="vegc"
typelist="recovery"

#explist="sba0056"
#vars="pft_nind pft_vegc"
#typelist="timmean AC${res}L05 AC${res}L01"


#typelist=""


## AC
#vars="vegc"
#typelist="AC${res}L01 SD${res}" # AC${res}L05 AC${res}L10 AC${res}L30"


## mesosouthamerica:
#explist="sba0080"
#vars="vegc agb fpc"
#typelist="timmean" # AC${res}L01 AC${res}L10"

# Nian 24 hypothese
#explist="sba0046"
#vars="agb_tree"
#typelist="timmean AC${res}L10" #AC${res}L01 

## sba0082-87
#explist="sba0082 sba0084 sba0087"
#explist="sba0080"
#vars="agb_tree agb vegc pft_nind pft_n_est fpc pft_vegc pft_chawo pft_cleaf pft_croot pft_csapw pft_bm_delta pft_wscal pft_bm_inc_carbon_afterallo pft_bm_inc_carbon_beforeallo pft_nind_beforemort pft_nind_beforeest pft_nind_afterest pft_est pft_mort pft_tinc_ind_cleaf pft_tinc_ind_croot pft_tinc_ind_csapw pft_tinc_ind_chawo"
#vars="pft_mort"
#typelist="timmean"

## runs with input2
##explist="sba0056 sba0061 sba0062 sba0063 sba0064 sba0065 sba0066 sba0067 sba0068"
##explist="sba0042 sba0044"
##vars="pft_nind_beforemort pft_bm_inc_carbon_beforeallo pft_n_est" 
##vars="agb pft_wscal swc1 swc2 swc pft_vegc pft_chawo pft_cleaf pft_croot pft_csapw pft_npp fpc"
#vars="pft_wscal"
#typelist="timmean"

#vars="lai_eff gpp npp fpc_woody" ## to check if resilience also works for other variables
#typelist="timmean AC${res}L05 AC${res}L01"

#vars="vegc agb cleafperm2 crootperm2 csapwperm2 chawoperm2 VOD1 VOD2 VOD3 VOD4 VOD5 VOD6"  # to check robustness of pseudo VOD
#typelist="timmean AC${res}L05 AC${res}L10"



#vars="vegc agb_tree agb"
#typelist="AC${res}L01 AC${res}L05 AC${res}L10 AC${res}L30 SD${res}"

## check agb_tree in other runs: => this is not the reason for the shifted MAPcrit
#explist="sba0042"
#vars="agb_tree"
#typelist="timmean AC${res}L10" #AC${res}L01 


###### AC of inputs
#explist="sba0044 sba0047"

## first: do this
#vars="pft_bm_inc_carbon_beforeallo pft_nind pft_wscal pft_nind_beforemort"
#typelist="timmean"

## second: run combine_variables to get pft_bm_inc_ind

## third: run this
#vars="pft_bm_inc_ind pft_nind pft_wscal"
#typelist="timmean AC${res}L01 AC${res}L10 AC${res}L30"

# sba0042: all PFTs
# sba0044: PFT1
# sba0047: PFT2
# sba0048: PFT3

## AC of vegc, and vegc per ind, and nind; each in total and per PFT:
#explist="sba0042 sba0044 sba0047 sba0048 sba0052 sba0053 sba0054"
##vars="pft_vegc pft_nind vegc nind vegcperind pft_vegcperind"
##vars="vegc nind agb"
#vars="vegc"
##typelist="timmean AC${res}L01 AC${res}L05 AC${res}L10 AC${res}L30"
#typelist="timmean timstd"

# est and mort budget
#vars="pft_est est pft_nind_beforemort pft_nind_aftermort pft_nind_beforeadj pft_nind_afteradj fpc"
#typelist="timmean"

#vars="pft_chawo pft_csapw pft_cleaf pft_croot pft_fracSH"
#typelist="timmean"



## diff exp
#explist="sba0053_minus_sba0044 sba0053_minus_sba0048 sba0054_minus_sba0047 sba0054_minus_sba0048"
#vars="est mort adjN adjNrel"
#typelist="timmean"


##############

#vars="pft_n_est"
#typelist="timmax mask"


## for calibration of the reduced model

#### before running this, run combine_variables.ksh
#vars="pft_P npp vegc agb litc pft_cleaf pft_csapw pft_croot pft_chawo fpc pft_vegc chawoperm2PLUSlit pft_chawoPLUSpft_sapw pft_chawoperm2PLUSpft_sapwperm2"

#vars="pft_nind"
#typelist="AC${res}L01 AC${res}L05 AC${res}L10"


#vars="litc cleaf csapw croot chawo agb vegc soilc chawoPLUSclit chawoPLUScsapw"
#typelist="timmean SD${res} AC${res}L01 AC${res}L05 AC${res}L10"

####


for var in ${vars}; do

  for exp in ${explist}; do

    mkdir -p ${exp}
    cd ${exp}
 
    if [[ ! -f ${var}.nc ]]; then
      echo "${var} not found"
      exit      
    fi

    for type in ${typelist}; do
    
     type_base2=`echo ${type} | cut -c1-2`
     
     if [[ ! -f ${var}_${type}.nc || ${calc_new} == 1 ]]; then      


      if [[ ${type} == "mask" ]]; then
          cdo timmax ${var}.nc ${var}_${type}_temp.nc
          cdo setrtoc,0.1,9e99,4 ${var}_${type}_temp.nc ${var}_${type}.nc
          rm ${var}_${type}_temp.nc
          
      # compute percentage of negative values
      elif [[ ${type} == yearnegcount ]]; then
          cdo yearsum ${var}.nc ${var}_yearsum.nc 
          cdo ltc,0 ${var}_yearsum.nc ${var}_lt0.nc 
          cdo timsum ${var}_lt0.nc ${var}_lt0sum.nc
          nofyears=`cdo showyear ${var}_yearsum.nc | wc -w` 
          cdo divc,${nofyears} -mulc,100 ${var}_lt0sum.nc ${var}_${type}.nc
          rm ${var}_lt0.nc ${var}_lt0sum.nc ${var}_yearsum.nc

      elif [[ ${type} == minmonmean ]]; then
        cdo yearmin ${var}.nc ${var}_yearmin.nc
        cdo timmean ${var}_yearmin.nc ${var}_${type}.nc
        rm ${var}_yearmin.nc

      elif [[ ${type} == freq0 ]]; then
        
        # frequency of zero increments
        nofyears=`cdo showyear D${var}.nc | wc -w`
        
        cdo gtc,1e-30  D${var}.nc D${var}_large.nc
        cdo ltc,-1e-30 D${var}.nc D${var}_small.nc 
        cdo add D${var}_large.nc D${var}_small.nc D${var}_one.nc

        cdo timsum D${var}_one.nc D${var}_sumofones.nc
        cdo divc,${nofyears} D${var}_sumofones.nc D${var}_freqC.nc
        cdo mulc,-1 -subc,1 D${var}_freqC.nc ${var}_${type}.nc

        rm -f D${var}_large.nc D${var}_small.nc
        rm -f D${var}_sumofones.nc D${var}_freqC.nc D${var}_one.nc
        cdo info ${var}_${type}.nc
        
      elif [[ ${type_base2} == "AC" || ${type_base2} == "SD" ]]; then   
        ### compute anomalies via ymonmean 
        if [[ ! -f ${var}_anomalies_${res}.nc ]]; then
     
          # remove mean (and annual cycle)
          if [[ ${res} == m ]]; then
          
            if [[ ! -f ${var}_ymonmean.nc ]]; then
              cdo ymonmean ${var}.nc ${var}_ymonmean.nc
            fi
            cdo sub ${var}.nc ${var}_ymonmean.nc ${var}_anomalies_${res}.nc
          else
            cdo sub ${var}.nc ${var}_timmean.nc ${var}_anomalies_${res}.nc
          fi #res
       
          # remove linear trend
          cdo detrend ${var}_anomalies_${res}.nc ${var}_anomalies_${res}.nc_detrended
          mv ${var}_anomalies_${res}.nc_detrended ${var}_anomalies_${res}.nc
       
        fi #anom
    
    
        if [[ ${type_base2} == SD ]]; then
          cdo timstd ${var}_anomalies_${res}.nc ${var}_${type}.nc

          # if 0, set to missing
          cdo setctomiss,0 ${var}_${type}.nc ${var}_${type}_miss.nc
          mv ${var}_${type}_miss.nc ${var}_${type}.nc        
    
        elif [[ ${type_base2} == AC ]]; then
          noftimesteps=`cdo showtime ${var}_anomalies_${res}.nc | wc -w` 
          
          lagall=`echo ${type} | cut -c5-6`
          lag1st=`echo ${type} | cut -c5`
          if [[ ${lag1st} == 0 ]]; then
            lag=`echo ${type} | cut -c6`
          else
            lag=${lagall}
          fi
    
          (( noftimesteps2 = ${noftimesteps} - ${lag} ))
          (( inisteps = 1 + ${lag} ))
        
          if [[ ${res} == m ]]; then
            cdo shifttime,${lag}month ${var}_anomalies_${res}.nc ${var}_shift.nc
          else
            cdo shifttime,${lag}year ${var}_anomalies_${res}.nc ${var}_shift.nc 
          fi
          cdo seltimestep,1/${noftimesteps2} ${var}_shift.nc ${var}_shift1.nc
          cdo seltimestep,${inisteps}/${noftimesteps} ${var}_anomalies_${res}.nc ${var}_shift2.nc 
          cdo timcor ${var}_shift1.nc ${var}_shift2.nc ${var}_${type}.nc
          
          # if -1, set to missing
          cdo setctomiss,-1 ${var}_${type}.nc ${var}_${type}_miss.nc
          mv ${var}_${type}_miss.nc ${var}_${type}.nc
          
          rm ${var}_shift1.nc ${var}_shift2.nc ${var}_shift.nc 
        fi #AC
    
    
        elif [[ ${type} == "recovery" ]]; then

          #export PYTHONPATH=$PYTHONPATH:/home/bathiany/Projects/Vegetation_resilience_indicators
          yr_perturb=0
          exp6digit=`echo ${exp} | cut -c6`
          exp7digit=`echo ${exp} | cut -c7`
          if [[ ${exp} == "sba0070" ]]; then
            exp_stationary="sba0056"
          elif [[ ${exp6digit} == "7" ]]; then
            exp_stationary="sba006${exp7digit}"
          else
            echo "stationary sim not defined"
          fi
          
          # compute anomaly w.r.t. stationary run (time dependent)
          years=$(cdo -s showyear ${var}.nc)
          start_year=$(cdo -s showyear ${var}.nc | head -c 5) # only works if start year has 4 digits
          start_year=$(echo ${start_year} | sed 's/ //g' ) 
          end_year=`echo ${years} | awk '{print $NF}'`
          cdo selyear,${start_year}/${end_year} ../${exp_stationary}/${var}.nc ${exp_stationary}_${var}.nc
          cdo sub ${var}.nc ${exp_stationary}_${var}.nc ${var}_anom2stat.nc

          cdo sub ${var}.nc ${exp_stationary}_${var}.nc ${var}_anom2stat.nc
          rm ${exp_stationary}_${var}.nc

          datapath=`pwd`
          #datapath="${datapath}/"
          echo "${var} ${datapath} ${yr_perturb}"
          
          python3 /home/bathiany/Projects/Vegetation_resilience_indicators/empirical_recovery_SB.py ${var} ${datapath} ${yr_perturb}


        
        #if [[ ( ${type} == yearmean || ${type} == timmax || ${type} == timmin || ${type} == timmean || ${type} == yearstd || ${type} == timstd || ${type} == ymonmin ) && -f ${var}.nc ]]; then
        else
          cdo ${type} ${var}.nc ${var}_${type}.nc
      
      fi #type
     fi #exists
    done #type

    cd ..

  done #var
done #exp

exit
