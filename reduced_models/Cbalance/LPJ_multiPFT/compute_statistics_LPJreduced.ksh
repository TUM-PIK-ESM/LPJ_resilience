#!/bin/ksh

vars="vegc"
explist="test1"
typelist="timmean ACyL01 ACyL02 ACyL05 ACyL10 ACyL20 ACyL30"

for exp in ${explist}; do
 cd /home/bathiany/Projects/Vegetation_resilience_indicators/reduced_models/Cbalance/LPJ_multiPFT/${exp}
 for var in ${vars}; do
  for type in ${typelist}; do

    origfile=${var}.nc
    targetfile=${var}_${type}.nc
    type_base2=`echo ${type} | cut -c1-2`

    if [[ ${type_base2} == "AC" || ${type_base2} == "SD" ]]; then  # type

        # anomalies, needs ymonmean or yearmean
        if [[ ! -f ${var}_anomalies.nc ]]; then # anom
         
          cdo sub ${var}.nc ${var}_timmean.nc ${var}_anomalies.nc

          # remove linear trend (not needed when stationary)
          cdo detrend ${var}_anomalies.nc ${var}_anomalies.nc_detrended
          mv ${var}_anomalies.nc_detrended ${var}_anomalies.nc
        fi #anom exists

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

    else
          
      cdo ${type} ${origfile} ${targetfile}
          
    fi #type

  done #type
 done #var
done #exp


exit
