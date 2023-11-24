#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Plot kinematic positions                                        ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 2.2                                                         ##
##                                                                           ##
##  DATE   : Mar-05, 2022                                                    ##
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
  if [ $# -ne 2 -a $# -ne 5 ]; then
    echo "#usage  : ./plotkin.sh kin_filename png_filename [x_ref y_ref z_ref]"
    echo "#example: ./plotkin.sh kin_2021149_rov1 rov1_2021149"
    echo "#if no x_ref y_ref z_ref, default x_avg y_avg z_avg"
    exit
  fi
  
  kinflname=$1
  pngflname=$2
  if [ $# -eq 5 ]; then
    x_ref=$3
    y_ref=$4
    z_ref=$5
  fi
  
  if [ $# -eq 5 ]; then
    enustd=`xyz2enu ${kinflname} enutmp ${x_ref} ${y_ref} ${z_ref}`
  else
    enustd=`xyz2enu ${kinflname} enutmp`
  fi
  rmse=`echo ${enustd} | awk '{printf("%16.2f", $2*100)}'`
  rmsn=`echo ${enustd} | awk '{printf("%16.2f", $3*100)}'`
  rmsu=`echo ${enustd} | awk '{printf("%16.2f", $4*100)}'`
  
  # time span
  mjd0=$(cat enutmp | head -n 1 | awk '{print $1}')
  mjd1=$(cat enutmp | tail -n 1 | awk '{print $1}')
  stime=$(cat enutmp | head -n 1 | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2}')
  etime=$(cat enutmp | tail -n 1 | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2}')
  dtime=$[$etime-$stime]
  
  # max min value
  min_e=`awk 'BEGIN{min =  999999999}{if($3 < min) min = $3;}END{print min}' enutmp`
  max_e=`awk 'BEGIN{max = -999999999}{if($3 > max) max = $3;}END{print max}' enutmp`
  max_e=`echo $max_e | awk '{print $1*100+0.5}'`
  min_e=`echo $min_e | awk '{print $1*100}'`
  min_n=`awk 'BEGIN{min =  999999999}{if($4 < min) min = $4;}END{print min}' enutmp`
  max_n=`awk 'BEGIN{max = -999999999}{if($4 > max) max = $4;}END{print max}' enutmp`
  max_n=`echo $max_n | awk '{print $1*100+0.5}'`
  min_n=`echo $min_n | awk '{print $1*100}'`
  min_u=`awk 'BEGIN{min =  999999999}{if($5 < min) min = $5;}END{print min}' enutmp`
  max_u=`awk 'BEGIN{max = -999999999}{if($5 > max) max = $5;}END{print max}' enutmp`
  max_u=`echo $max_u | awk '{print $1*100+0.5}'`
  min_u=`echo $min_u | awk '{print $1*100}'`
  min_dop=`awk 'BEGIN{min =  999999999}{if($13 < min) min = $13;}END{print min}' enutmp`
  max_dop=`awk 'BEGIN{max = -999999999}{if($13 > max) max = $13;}END{print max}' enutmp`
  min_sat=`awk 'BEGIN{min =  999999999}{if($6 < min) min = $6;}END{print min}' enutmp`
  max_sat=`awk 'BEGIN{max = -999999999}{if($6 > max) max = $6;}END{print max}' enutmp`
  min_dop=`echo $min_dop | awk '{print $1-0.5}'`
  max_dop=`echo $max_dop | awk '{print $1+1.5}'`
  ((min_sat=min_sat-1))
  ((max_sat=max_sat+2))
  
  gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,black
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
  
  PS="${pngflname}.ps"
  J="X18c/3c"
  Re="${stime}/${etime}/${min_e}/${max_e}"
  Rn="${stime}/${etime}/${min_n}/${max_n}"
  Ru="${stime}/${etime}/${min_u}/${max_u}"
  Rdop="${stime}/${etime}/${min_dop}/${max_dop}"
  Rnum="${stime}/${etime}/${min_sat}/${max_sat}"

  # Time axis
  readonly fmt="%04d-%02d-%02dT%02d:%02d:%02d"
  sod0=$(cat enutmp | head -n 1 | awk '{print $2}')
  sod1=$(cat enutmp | tail -n 1 | awk '{print $2}')
  ymd0=$(ydoy2ymd $(mjd2ydoy $mjd0))
  ymd1=$(ydoy2ymd $(mjd2ydoy $mjd1))
  stime=$(echo $ymd0 $sod0 | awk -v fmt=$fmt '{hh=int($4/3600);mm=int($4/60)-hh*60;ss=$4%60;printf(fmt,$1,$2,$3,hh,mm,ss)}')
  etime=$(echo $ymd1 $sod1 | awk -v fmt=$fmt '{hh=int($4/3600);mm=int($4/60)-hh*60;ss=$4%60;printf(fmt,$1,$2,$3,hh,mm,ss)}')

  Rte="${stime}/${etime}/${min_e}/${max_e}"
  Rtn="${stime}/${etime}/${min_n}/${max_n}"
  Rtu="${stime}/${etime}/${min_u}/${max_u}"
  Rtd="${stime}/${etime}/${min_dop}/${max_dop}"

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

  gmt psxy -J$J -R${Ru} -T -K > $PS
  
  # E N U
  gmt psbasemap -J$J -R${Re} -BWe -Byaf+l"East (cm)" -Y11c -K -O >> $PS
  cat enutmp | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$3*100}' | gmt psxy -J$J -R${Re} -Fas -W0.03c -K -O >> $PS
  echo East  ${rmse} cm | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PS
  gmt psbasemap -J$J -R${Rte} -Bsn $Bx1 -K -O >> $PS
  
  gmt psbasemap -J$J -R${Rn} -BWe -Byaf+l"North (cm)" -Y-3.2c -K -O >> $PS
  cat enutmp | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$4*100}' | gmt psxy -J$J -R${Rn} -Fas -W0.03c -K -O >> $PS
  echo North  ${rmsn} cm | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PS
  gmt psbasemap -J$J -R${Rtn} -Bsn $Bx1 -K -O >> $PS
  
  gmt psbasemap -J$J -R${Ru} -BWe -Byaf+l"Up (cm)" -Y-3.2c -K -O >> $PS
  cat enutmp | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$5*100}' | gmt psxy -J$J -R${Ru} -Fas -W0.03c -K -O >> $PS
  echo Up  ${rmsu} cm | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PS
  gmt psbasemap -J$J -R${Rtu} -Bsn $Bx1 -K -O >> $PS
  
  # DOP & Satnum
  gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,red
  gmt set FONT_LABEL 12p,Helvetica,red
  gmt set MAP_FRAME_PEN thinner,red 
  gmt set MAP_TICK_PEN thinner,red
  gmt psbasemap -J$J -R${Rdop} -BW -Byaf+l"DOP" -Y-3.2c -K -O >> $PS
  gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,blue
  gmt set FONT_LABEL 12p,Helvetica,blue
  gmt set MAP_FRAME_PEN thinner,blue
  gmt set MAP_TICK_PEN thinner,blue
  gmt psbasemap -J$J -R${Rnum} -BE -Byaf+l"Satellite number" -K -O >> $PS
  cat enutmp | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$13}' | gmt psxy -J$J -R${Rdop} -Fas -W0.03c,red  -K -O >> $PS
  cat enutmp | awk -v mjd=$mjd0 '{print ($1-mjd)*86400+$2,$6}' | gmt psxy -J$J -R${Rnum} -Fas -W0.03c,blue -K -O >> $PS
  gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,black
  gmt set FONT_ANNOT_SECONDARY 11p,Helvetica,black
  gmt set FONT_LABEL 12p,Helvetica,black
  gmt set MAP_FRAME_PEN thinner,black
  gmt set MAP_TICK_PEN thinner,black
  gmt psbasemap -J$J -R${Rtd} -BSn $Bx1+l"Time" $Bx2 -K -O >> $PS

  > legend.txt
  echo "N 2" >> legend.txt
  echo "S -3.3c c 0.1c red  1p,red,solid  -3.1c PDOP"             >> legend.txt
  echo "S -2.5c c 0.1c blue 1p,blue,solid -2.3c Satellite number" >> legend.txt
  cat legend.txt | gmt pslegend -J -R -DjTR+w1.5c+l0.8+o0.1c/0.1c  -K -O >> $PS
  
  gmt psxy -J$J -R${Rnum} -T -O >> $PS
  
  # gmt clear history
  rm *.history
  gmt psconvert -Tg -A -P -Z $PS
  
  rm -f gmt.*  *.eps *.ps
  rm -f enutmp legend.txt
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
