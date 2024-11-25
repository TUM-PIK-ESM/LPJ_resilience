#!/bin/ksh

### takes basic output from reduced model and combines it to derive more variables

#exp="test1"  # 1, 5
#exp_driving="sba0056"
#NPFT=8

#exp="test2" # 2, 6
#exp_driving="sba0061"
#NPFT=1

#exp="test4"
#exp_driving="sba0042"
#NPFT=8

#exp="test9" # test 3, 7, 8, 9
#exp_driving="sba0044"
#NPFT=1

exp_list="sba0056_P0000 sba0056_P0010 sba0056_P0011 sba0056_P0012 sba0056_P0013 sba0056_P0014 sba0056_P0015 sba0056_P0020 sba0056_P0021 sba0056_P0022 sba0056_P0023 sba0056_P0024 sba0056_P0025 sba0056_P0030 sba0056_P0031 sba0056_P0032 sba0056_P0033 sba0056_P0034 sba0056_P0035"
exp_list="sba0056_P0041"

exp_driving="sba0056"
NPFT=8

calc_new=0


for exp in ${exp_list}; do

cd /home/bathiany/Projects/Vegetation_resilience_indicators/reduced_models/Cbalance/LPJ_multiPFT/${exp}

### carbon pools per m2: g/ind x ind/m2
for pool in cleaf csapw croot chawo; do
 if [[ ! -f ${pool}perm2.nc || ${calc_new} == 1 ]]; then
    cdo mul pft_nind.nc pft_${pool}.nc pft_${pool}perm2.nc
    cdo timmean pft_${pool}perm2.nc pft_${pool}perm2_timmean.nc
    cdo vertsum pft_${pool}perm2.nc ${pool}perm2.nc

    if [[ ${pool} == "chawo" ]]; then
        varname="C_heartwood_pft"
    elif [[ ${pool} == "csapw" ]]; then
        varname="C_sapwood_pft"
    elif [[ ${pool} == "croot" ]]; then
        varname="C_root_pft"
    elif [[ ${pool} == "cleaf" ]]; then
        varname="C_leaf_pft"
    fi

    cdo setname,${varname} ${pool}perm2.nc named.nc
    mv named.nc ${pool}perm2.nc

 fi
done

## vegc
if [[ ! -f pft_vegc.nc || ${calc_new} == 1 ]]; then
 cdo add pft_cleafperm2.nc pft_crootperm2.nc LR.nc
 cdo add pft_csapwperm2.nc pft_chawoperm2.nc SH.nc
 cdo add LR.nc SH.nc pft_vegc.nc
 rm LR.nc SH.nc
 cdo setname,VegC_pft pft_vegc.nc named.nc
 mv named.nc pft_vegc.nc
fi


if [[ ! -f vegc.nc || ${calc_new} == 1 ]]; then
  cdo vertsum pft_vegc.nc vegc.nc
  cdo setname,VegC vegc.nc named.nc
  mv named.nc vegc.nc
fi

## agb
if [[ ! -f pft_agb.nc || ${calc_new} == 1 ]]; then
 cdo mulc,0.66 pft_csapwperm2.nc sap66.nc
 cdo add pft_chawoperm2.nc pft_cleafperm2.nc pft_LH.nc
 cdo add pft_LH.nc sap66.nc pft_agb.nc
 rm sap66.nc pft_LH.nc
fi

if [[ ! -f agb.nc || ${calc_new} == 1 ]]; then
  cdo vertsum pft_agb.nc agb.nc
fi

if [[ ! -f nind.nc || ${calc_new} == 1 ]]; then
  cdo vertsum pft_nind.nc nind.nc
fi



done


exit
#### eval with LPJ
#varlist="pft_cleaf pft_croot pft_csapw pft_chawo pft_vegc pft_nind fpc vegc"
varlist="vegc agb"

LPJdir="/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/"

for exp in ${exp_list}; do
for var in ${varlist}; do
  
  years=$(cdo -s showyear ${var}.nc)
  start_year=$(cdo -s showyear ${var}.nc | head -c 5) # only works if start year has 4 digits
  start_year=$(echo ${start_year} | sed 's/ //g' ) 
  end_year=`echo ${years} | awk '{print $NF}'`
  var_first3=`echo ${var} | head -c 3`

  if [[ ! -f  ${var}_diff2LPJ_rmse.nc || ${calc_new} == 1 ]]; then

   if [[ ${var_first3} == "pft" || ${var_first3} == "fpc" ]]; then
      cdo -s sellevel,1/${NPFT} -selyear,${start_year}/${end_year} ${LPJdir}/${exp_driving}/${var}.nc \
                           ${var}_LPJ.nc
    else
      cdo -s selyear,${start_year}/${end_year} ${LPJdir}/${exp_driving}/${var}.nc \
                           ${var}_LPJ.nc
   fi

   cdo -s timcor ${var}.nc ${var}_LPJ.nc ${var}_corr_LPJ.nc

   cdo -s sub ${var}.nc ${var}_LPJ.nc ${var}_diff2LPJ.nc
   cdo timmean ${var}_diff2LPJ.nc ${var}_bias_LPJ.nc

   cdo sqr ${var}_diff2LPJ.nc ${var}_diff2LPJ_sqr.nc
   cdo timmean ${var}_diff2LPJ_sqr.nc ${var}_diff2LPJ_sqr_sum.nc
   cdo sqrt ${var}_diff2LPJ_sqr_sum.nc ${var}_rmse_LPJ.nc
   rm ${var}_diff2LPJ_sqr_sum.nc ${var}_diff2LPJ_sqr.nc

   ## use est mask to remove all types that are not established:
   if [[ ${var_first3} == "pft" || ${var_first3} == "fpc" ]]; then
     cdo -s sellevel,1/${NPFT} ${LPJdir}/${exp_driving}/pft_est_timmax.nc mask.nc
   else
     cdo -s sellevel,1 ${LPJdir}/${exp_driving}/pft_est_timmax.nc mask.nc
   fi
   
   cdo ifthen mask.nc ${var}_corr_LPJ.nc ${var}_corr_LPJ.nc_2
   mv ${var}_corr_LPJ.nc_2 ${var}_corr_LPJ.nc
   rm mask.nc
   
  fi
  
done #var

done #exp


exit

