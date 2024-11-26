#!/bin/ksh

exp_list="sba0056_P0000"

#vars="vegc pft_vegc pft_cleafperm2 cleafperm2 pft_crootperm2 crootperm2 pft_csapwperm2 csapwperm2 pft_chawoperm2 chawoperm2"
vars="vegc pft_vegc nind pft_nind"
typelist="timmean ACyL01 ACyL02 ACyL05 ACyL10 ACyL20 ACyL30"

#vars="vegc"
#typelist="timmean ACyL05"

exp_list="sba0056_P0010 sba0056_P0011 sba0056_P0012 sba0056_P0013 sba0056_P0014 sba0056_P0015 sba0056_P0020 sba0056_P0021 sba0056_P0022 sba0056_P0023 sba0056_P0024 sba0056_P0025 sba0056_P0030 sba0056_P0031 sba0056_P0032 sba0056_P0033 sba0056_P0034 sba0056_P0035 sba0056_P0040 sba0056_P0041 sba0056_P0042 sba0056_P0043 sba0056_P0044"

exp_list="sba0056_P0010"

vars="vegc"
typelist="recovery"
yr_perturb=200


for exp in ${exp_list}; do
 cd /home/bathiany/Projects/Vegetation_resilience_indicators/reduced_models/Cbalance/LPJ_multiPFT/${exp}
 exp_stationary_head=`echo ${exp} | cut -c1-7`
 exp_stationary="${exp_stationary_head}_P0000"
 for var in ${vars}; do
  for type in ${typelist}; do

    origfile=${var}.nc
    targetfile=${var}_${type}.nc
    type_base2=`echo ${type} | cut -c1-2`

    if [[ ${type_base2} == "AC" || ${type_base2} == "SD" ]]; then  # type

        # anomalies, needs ymonmean or yearmean
        #if [[ ! -f ${var}_anomalies.nc ]]; then # anom
         
          cdo sub ${var}.nc ${var}_timmean.nc ${var}_anomalies.nc

          # remove linear trend (not needed when stationary)
          cdo detrend ${var}_anomalies.nc ${var}_anomalies.nc_detrended
          mv ${var}_anomalies.nc_detrended ${var}_anomalies.nc
        #fi #anom exists

        # stdev, needs anomalies
        if [[ ${type_base2} == SD ]]; then
          cdo timstd ${var}_anomalies.nc ${targetfile}

          # if timstd 0, set to missing
          cdo setctomiss,0 ${targetfile} ${var}_${type}_miss.nc
          mv ${var}_${type}_miss.nc ${targetfile}        

        # AC, needs anomalies
        elif [[ ${type_base2} == AC ]]; then
          noftimesteps=`cdo showtime ${var}_anomalies.nc | wc -w` 
          
          lagall=`echo ${type} | cut -c5-6`
          lag1st=`echo ${type} | cut -c5`
          if [[ ${lag1st} == 0 ]]; then
            lag=`echo ${type} | cut -c6`
          else
            lag=${lagall}
          fi
    
          (( noftimesteps2 = ${noftimesteps} - ${lag} ))
          (( inisteps = 1 + ${lag} ))
        
          cdo shifttime,${lag}year ${var}_anomalies.nc ${var}_shift.nc 
          
          cdo seltimestep,1/${noftimesteps2} ${var}_shift.nc ${var}_shift1.nc
          cdo seltimestep,${inisteps}/${noftimesteps} ${var}_anomalies.nc \
                                                      ${var}_shift2.nc 
          cdo timcor ${var}_shift1.nc ${var}_shift2.nc ${targetfile}
          
          # if -1, set to missing
          cdo setctomiss,-1 ${targetfile} ${var}_${type}_miss.nc
          mv ${var}_${type}_miss.nc ${targetfile}
          
          rm ${var}_shift1.nc ${var}_shift2.nc ${var}_shift.nc 
        fi #AC

    elif [[ ${type} == "recovery" ]]; then  # type

     if [[ ! -f ${var}_recoveryrate.nc ]]; then 
      ## use recovery_template code (old)
      #cp ../recovery_template.py recovery_${var}_${exp}.py
            
      ## compute anomaly w.r.t. stationary run's mean
      #cdo timmean ../${exp_stationary}/${var}.nc ${exp_stationary}_${var}.nc
      
      # compute anomaly w.r.t. stationary run (time dependent)
      years=$(cdo -s showyear ${var}.nc)
      start_year=$(cdo -s showyear ${var}.nc | head -c 5) # only works if start year has 4 digits
      start_year=$(echo ${start_year} | sed 's/ //g' ) 
      end_year=`echo ${years} | awk '{print $NF}'`
      cdo selyear,${start_year}/${end_year} ../${exp_stationary}/${var}.nc ${exp_stationary}_${var}.nc
      cdo sub ${var}.nc ${exp_stationary}_${var}.nc ${var}_anom2stat.nc

      cdo sub ${var}.nc ${exp_stationary}_${var}.nc ${var}_anom2stat.nc
      rm ${exp_stationary}_${var}.nc

      ## modify content: (old)
      #sed -i "s/experiment_id/${exp}/g" recovery_${var}_${exp}.py
      #sed -i "s/var_id/${var}/g" recovery_${var}_${exp}.py
      #sed -i "s/yr_perturb_id/${yr_perturb}/g" recovery_${var}_${exp}.py
      ##exit
      #python3 recovery_${var}_${exp}.py
      #rm recovery_${var}_${exp}.py   
      datapath=`pwd`
      python3 /home/bathiany/Projects/Vegetation_resilience_indicators/empirical_recovery_SB.py ${var} ${datapath} ${yr_perturb}

     fi
    else
          
      cdo ${type} ${origfile} ${targetfile}
          
    fi #type

  done #type
 done #var
done #exp


exit
