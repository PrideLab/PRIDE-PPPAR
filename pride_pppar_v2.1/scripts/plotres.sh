#!/bin/bash
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

gmt set FONT_ANNOT_PRIMARY 11p,Helvetica,black
gmt set FONT_LABEL 12p,Helvetica,black
gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
gmt set MAP_LABEL_OFFSET 0.1c
gmt set MAP_TICK_LENGTH_PRIMARY -0.08c
gmt set MAP_FRAME_WIDTH 0.001c
gmt set MAP_FRAME_PEN thinner,black
gmt set MAP_GRID_PEN_PRIMARY 0.2p,darkgray,-
gmt set FONT_TITLE 20p,Helvetica,black
gmt set MAP_TITLE_OFFSET 0.1c

#time span
mjd0=`cat $res_flname | head -n 1 | awk '{print $1}'`
stime=`cat $res_flname | head -n 1 | awk '{print ($1-$mjd0)*86400+$2}'`
etime=`cat $res_flname | tail -n 1 | awk '{print ($1-$mjdo)*86400+$2}'`
#max min value
min_p=`awk 'BEGIN{min = 999999999}{if($4 < min) min = $4;}END{print min}' $res_flname`
max_p=`awk 'BEGIN{max = -999999999}{if($4 > max) max = $4;}END{print max}' $res_flname`
min_r=`awk 'BEGIN{min = 999999999}{if($5 < min) min = $5;}END{print min}' $res_flname`
max_r=`awk 'BEGIN{max = -999999999}{if($5 > max) max = $5;}END{print max}' $res_flname`
echo $min_p $max_p $min_r $max_r

PS="${pngflname}.ps"
J=X18c/5c
Rp=${stime}/${etime}/${min_p}/${max_p}
Rr=${stime}/${etime}/${min_r}/${max_r}

#all sat's res in one picture 
gmt psxy -J$J -R${Rp} -T -K > $PS

#Phase
#gmt psbasemap -J$J -R${Rp} -BWSen -Bxaf -Byaf+l"Carrier-Phase Residuals(m)" -Y11c -K -O >> $PS
gmt psbasemap -J$J -R${Rp} -BWsen -Bxaf -Byaf+l"Carrier-Phase Residuals(m)" -Y10.5c -K -O >> $PS
for prn in ${prns[*]}
do
  lnum=`cat $res_flname | grep "${prn}" | wc -l`
  [ ${lnum} -eq 0 ] && continue
  cat $res_flname | grep "${prn}" | awk '{print ($1-$mjd0)*86400+$2,$4}' | gmt psxy -J$J -R${Rp} -Sc0.05c -Gblack -K -O >> $PS
done

#Range
gmt psbasemap -J$J -R${Rr} -BWSen -Bxaf+l"Seconds of day (s)" -Byaf+l"Pseudorange Residuals(m)" -Y-5.2c -K -O >> $PS
for prn in ${prns[*]}
do
  lnum=`cat $res_flname | grep "${prn}" | wc -l`
  [ ${lnum} -eq 0 ] && continue
  cat $res_flname | grep "${prn}" | awk '{print ($1-$mjd0)*86400+$2,$5}' | gmt psxy -J$J -R${Rr} -Sc0.05c -Gblack -K -O >> $PS
done

gmt psxy -J$J -R${Rr} -T -O >> $PS

#Single Satellite
for prn in ${prns[*]}
do
  lnum=`cat $res_flname | grep "${prn}" | wc -l`
  [ ${lnum} -eq 0 ] && continue
  echo single satellite:$prn 

  PSs="${pngflname}_${prn}.ps"
  gmt psxy -J$J -R${Rp} -T -K > $PSs
  
  gmt psbasemap -J$J -R${Rp} -BWsen -Bxaf -Byaf+l"Carrier-Phase Residuals(m)" -Y10.5c -K -O >> $PSs
  cat $res_flname | grep "${prn}" | awk '{print ($1-$mjd0)*86400+$2,$4}' | gmt psxy -J$J -R${Rp} -Sc0.05c -Gblack -K -O >> $PSs
  echo $prn | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PSs

  gmt psbasemap -J$J -R${Rr} -BWSen -Bxaf+l"Seconds of day (s)" -Byaf+l"Pseudorange Residuals(m)" -Y-5.2c -K -O >> $PSs
  cat $res_flname | grep "${prn}" | awk '{print ($1-$mjd0)*86400+$2,$5}' | gmt psxy -J$J -R${Rr} -Sc0.05c -Gblack -K -O >> $PSs
  
  gmt psxy -J$J -R${Rr} -T -O >> $PSs
  gmt psconvert -Tg -A -P -Z  $PSs
done


#gmt clear history
rm *.history
echo all satellites
gmt psconvert -Tg -A -P -Z $PS

rm -f gmt.*  *.eps *.ps


