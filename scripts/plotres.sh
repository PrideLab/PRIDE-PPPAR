#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Plot pseudorange and carrier-phase residuals                    ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 2.2                                                         ##
##                                                                           ##
##  DATE   : Jan-18, 2022                                                    ##
##                                                                           ##
##              @ GNSS RESEARCH CENTER, WUHAN UNIVERSITY, 2023               ##
##                                                                           ##
##    Copyright (C) 2022 by Wuhan University                                 ##
##                                                                           ##
##    This program is free software: you can redistribute it and/or modify   ##
##    it under the terms of the GNU General Public License (version 3) as    ##
##    published by the Free Software Foundation.                             ##
##                                                                           ##
##    This program is distributed in the hope that it will be useful,        ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of         ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           ##
##    GNU General Public License (version 3) for more details.               ##
##                                                                           ##
##    You should have received a copy of the GNU General Public License      ##
##    along with this program. If not, see <https://www.gnu.org/licenses/>.  ##
##                                                                           ##
###############################################################################

main() {
  if [ $# -ne 2 ]; then
    echo "#usage  : ./plotres.sh res_filename png_filename "
    echo "#example: ./plotres.sh res_2021149_rov1 rov1_res_2021149 "
    exit
  fi
  
  resflname=$1
  pngflname=$2
  res_flname="${resflname}_"
  > $res_flname 
  cat $resflname | while read line
  do
    if [[ ${line:0:3} == "TIM" ]];then
      mjd=`echo $line | awk '{print $8}'`
      sod=`echo $line | awk '{print $9}'`
      continue
    fi
    echo  $mjd $sod $line >> $res_flname
  done
  nhead=`grep -n "END OF HEADER" $res_flname | cut -d ":" -f 1`
  ((nhead=nhead+1))
  sed -i "1,${nhead}d" $res_flname
  
  #plot SAT[..] IONPS
  prns=(G01 G02 G03 G04 G05 G06 G07 G08 G09 G10 G11 G12 G13 G14 G15 G16 \
        G17 G18 G19 G20 G21 G22 G23 G24 G25 G26 G27 G28 G29 G30 G31 G32 \
        R01 R02 R03 R04 R05 R06 R07 R08 R09 R10 R11 R12 R13 R14 R15 R16 \
        R17 R18 R19 R20 R21 R22 R23 R24 E01 E02 E03 E04 E05 E06 E07 E08 \
        E09 E10 E11 E12 E13 E14 E15 E16 E17 E18 E19 E20 E21 E22 E23 E24 \
        E25 E26 E27 E28 E29 E30 E31 E32 E33 E34 E35 E36 C01 C02 C03 C04 \
        C05 C06 C07 C08 C09 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 \
        C21 C22 C23 C24 C25 C26 C27 C28 C29 C30 C31 C32 C33 C34 C35 C36 \
        C37 C38 C39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C56 C57 C58 C59 \
        C60 C61 J01 J02 J03 J07)
  
  gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,black
  gmt set FONT_ANNOT_SECONDARY 11p,Helvetica,black
  gmt set FONT_LABEL 12p,Helvetica,black
  gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
  gmt set MAP_LABEL_OFFSET 0.1c
  gmt set MAP_TICK_LENGTH_PRIMARY -0.08c
  gmt set MAP_FRAME_WIDTH 0.001c
  gmt set MAP_FRAME_PEN thinner,black
  gmt set MAP_GRID_PEN_PRIMARY 0.2p,darkgray,-
  gmt set FONT_TITLE 20p,Helvetica,black
  gmt set MAP_TITLE_OFFSET 0.1c
  gmt set FORMAT_DATE_MAP yyyy-mm-dd FORMAT_CLOCK_MAP hh:mm
  
  # time span
  mjd0=$(cat $res_flname | head -n 1 | awk '{print $1}')
  mjd1=$(cat $res_flname | tail -n 1 | awk '{print $1}')
  stime=$(cat $res_flname | head -n 1 | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2}')
  etime=$(cat $res_flname | tail -n 1 | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2}')
  dtime=$[$etime-$stime]
  
  # max min value
  min_p=`awk 'BEGIN{min =  999999999}{if($4 < min) min = $4;}END{print min}' $res_flname`
  max_p=`awk 'BEGIN{max = -999999999}{if($4 > max) max = $4;}END{print max}' $res_flname`
  min_r=`awk 'BEGIN{min =  999999999}{if($5 < min) min = $5;}END{print min}' $res_flname`
  max_r=`awk 'BEGIN{max = -999999999}{if($5 > max) max = $5;}END{print max}' $res_flname`
  echo $min_p $max_p $min_r $max_r
  
  PS="${pngflname}.ps"
  J="X18c/5c"
  Rp="${stime}/${etime}/${min_p}/${max_p}"
  Rr="${stime}/${etime}/${min_r}/${max_r}"

  # Time axis
  readonly fmt="%04d-%02d-%02dT%02d:%02d:%02d"
  sod0=$(cat $res_flname | head -n 1 | awk '{print $2}')
  sod1=$(cat $res_flname | tail -n 1 | awk '{print $2}')
  ymd0=$(ydoy2ymd $(mjd2ydoy $mjd0))
  ymd1=$(ydoy2ymd $(mjd2ydoy $mjd1))
  stime=$(echo $ymd0 $sod0 | awk -v fmt=$fmt '{hh=int($4/3600);mm=int($4/60)-hh*60;ss=$4%60;printf(fmt,$1,$2,$3,hh,mm,ss)}')
  etime=$(echo $ymd1 $sod1 | awk -v fmt=$fmt '{hh=int($4/3600);mm=int($4/60)-hh*60;ss=$4%60;printf(fmt,$1,$2,$3,hh,mm,ss)}')

  Rtp="${stime}/${etime}/${min_p}/${max_p}"
  Rtr="${stime}/${etime}/${min_r}/${max_r}"

  if [ $dtime -le 86400 ]; then
    mp=$[$[$[dtime+10799]/10800]*15]
    ms=$[$[$[dtime+10799]/10800]*3]
    Bx1="-Bpxa${mp}Mf${ms}m"
    if [ $mp -ge 60 ]; then
      hp=$[$[$mp+59]/60]
      ms=30
      Bx1="-Bpxa${hp}Hf${ms}m"
    fi
  else
    hp=$[12/$[$mjd1-$mjd0+1]]
    hp=$[$[$hp+23]/$hp]
    hs=1
    Bx1="-Bpxa${hp}Hf${hs}h"
  fi

  [ $mjd0 -ne $mjd1 ] && Bx2="-Bsxa1D"
  
  # all sat's res in one picture 
  gmt psxy -J$J -R${Rp} -T -K > $PS
  
  # Phase
  gmt psbasemap -J$J -R${Rp} -BWe -Byaf+l"Carrier-Phase Residuals(m)" -Y10.5c -K -O >> $PS
  for prn in ${prns[*]}
  do
    lnum=`cat $res_flname | grep "${prn}" | wc -l`
    [ ${lnum} -eq 0 ] && continue
    cat $res_flname | grep "${prn}" | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$4}' | gmt psxy -J$J -R${Rp} -Sc0.05c -Gblack -K -O >> $PS
  done
  gmt psbasemap -J$J -R${Rtp} -Bsn $Bx1 -K -O >> $PS
  
  # Range
  gmt psbasemap -J$J -R${Rr} -BWe -Byaf+l"Pseudorange Residuals(m)" -Y-5.2c -K -O >> $PS
  for prn in ${prns[*]}
  do
    lnum=`cat $res_flname | grep "${prn}" | wc -l`
    [ ${lnum} -eq 0 ] && continue
    cat $res_flname | grep "${prn}" | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$5}' | gmt psxy -J$J -R${Rr} -Sc0.05c -Gblack -K -O >> $PS
  done
  gmt psbasemap -J$J -R${Rtr} -BSn $Bx1+l"Time" -K -O >> $PS
  
  gmt psxy -J$J -R${Rr} -T -O >> $PS
  
  # Single Satellite
  for prn in ${prns[*]}
  do
    lnum=`cat $res_flname | grep "${prn}" | wc -l`
    [ ${lnum} -eq 0 ] && continue
    echo "single satellite: $prn"
  
    PSs="${pngflname}_${prn}.ps"
    gmt psxy -J$J -R${Rp} -T -K > $PSs
    
    gmt psbasemap -J$J -R${Rp} -BWe -Byaf+l"Carrier-Phase Residuals(m)" -Y10.5c -K -O >> $PSs
    cat $res_flname | grep "${prn}" | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$4}' | gmt psxy -J$J -R${Rp} -Sc0.05c -Gblack -K -O >> $PSs
    echo $prn | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PSs
    gmt psbasemap -J$J -R${Rtp} -Bsn $Bx1 -K -O >> $PSs
  
    gmt psbasemap -J$J -R${Rr} -BWe -Byaf+l"Pseudorange Residuals(m)" -Y-5.2c -K -O >> $PSs
    cat $res_flname | grep "${prn}" | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$5}' | gmt psxy -J$J -R${Rr} -Sc0.05c -Gblack -K -O >> $PSs
    gmt psbasemap -J$J -R${Rtr} -BSn $Bx1+l"Time" $Bx2 -K -O >> $PSs
    
    gmt psxy -J$J -R${Rr} -T -O >> $PSs
    gmt psconvert -Tg -A -P -Z  $PSs
  done
  
  # gmt clear history
  rm *.history
  echo "all satellites"
  gmt psconvert -Tg -A -P -Z $PS
  
  rm -f gmt.*  *.eps *.ps
}

######################################################################
##                      Time Convert Funcitons                      ##
##  Author: Shuyin Mao      shuyinm@whu.edu.cn                      ##
######################################################################
ymd2mjd() {
    local year=$1
    local mon=$((10#$2))
    local day=$((10#$3))
    [ $year -lt 100 ] && year=$((year+2000))
    if [ $mon -le 2 ];then
        mon=$(($mon+12))
        year=$(($year-1))
    fi
    local mjd=`echo $year | awk '{print $1*365.25-$1*365.25%1-679006}'`
    mjd=`echo $mjd $year $mon $day | awk '{print $1+int(30.6001*($3+1))+2-int($2/100)+int($2/400)+$4}'`
    echo $mjd
}

mjd2ydoy() {
    local mjd=$1
    local year=$((($mjd + 678940)/365))
    local mjd0=$(ymd2mjd $year 1 1)
    local doy=$(($mjd-$mjd0))
    while [ $doy -le 0 ];do
        year=$(($year-1))
        mjd0=$(ymd2mjd $year 1 1)
        doy=$(($mjd-$mjd0+1))
    done
    printf "%d %03d\n" $year $doy
}

ydoy2ymd() {
    local iyear=$1
    local idoy=$((10#$2))
    local days_in_month=(31 28 31 30 31 30 31 31 30 31 30 31)
    local iday=0
    [ $iyear -lt 100 ] && iyear=$((iyear+2000))
    local tmp1=$(($iyear%4))
    local tmp2=$(($iyear%100))
    local tmp3=$(($iyear%400))
    if [ $tmp1 -eq 0 -a $tmp2 -ne 0 ] || [ $tmp3 -eq 0 ]; then
       days_in_month[1]=29
    fi
    local id=$idoy
    local imon=0
    local days
    for days in ${days_in_month[*]}
    do
        id=$(($id-$days))
        imon=$(($imon+1))
        if [ $id -gt 0 ]; then
            continue
        fi
        iday=$(($id + $days))
        break
    done
    printf "%d %02d %02d\n" $iyear $imon $iday
}

main "$@"

