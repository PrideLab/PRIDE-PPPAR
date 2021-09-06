#!/bin/bash

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


#time span
mjd0=`cat enutmp | head -n 1 | awk '{print $1}'`
stime=`cat enutmp | head -n 1 | awk '{print ($1-$mjd0)*86400+$2}'`
etime=`cat enutmp | tail -n 1 | awk '{print ($1-$mjd0)*86400+$2}'`

#max min value
min_e=`awk 'BEGIN{min = 999999999}{if($3 < min) min = $3;}END{print min}' enutmp`
max_e=`awk 'BEGIN{max = -999999999}{if($3 > max) max = $3;}END{print max}' enutmp`
max_e=`echo $max_e | awk '{print $1*100+0.5}'`
min_e=`echo $min_e | awk '{print $1*100}'`
min_n=`awk 'BEGIN{min = 999999999}{if($4 < min) min = $4;}END{print min}' enutmp`
max_n=`awk 'BEGIN{max = -999999999}{if($4 > max) max = $4;}END{print max}' enutmp`
max_n=`echo $max_n | awk '{print $1*100+0.5}'`
min_n=`echo $min_n | awk '{print $1*100}'`
min_u=`awk 'BEGIN{min = 999999999}{if($5 < min) min = $5;}END{print min}' enutmp`
max_u=`awk 'BEGIN{max = -999999999}{if($5 > max) max = $5;}END{print max}' enutmp`
max_u=`echo $max_u | awk '{print $1*100+0.5}'`
min_u=`echo $min_u | awk '{print $1*100}'`
min_dop=`awk 'BEGIN{min = 999999999}{if($7 < min) min = $7;}END{print min}' enutmp`
max_dop=`awk 'BEGIN{max = -999999999}{if($7 > max) max = $7;}END{print max}' enutmp`
min_sat=`awk 'BEGIN{min = 999999999}{if($6 < min) min = $6;}END{print min}' enutmp`
max_sat=`awk 'BEGIN{max = -999999999}{if($6 > max) max = $6;}END{print max}' enutmp`
min_dop=`echo $min_dop | awk '{print $1-0.5}'`
max_dop=`echo $max_dop | awk '{print $1+1.5}'`
((min_sat=min_sat-1))
((max_sat=max_sat+2))


gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,black
gmt set FONT_LABEL 11p,Helvetica,black
gmt set MAP_ANNOT_OFFSET_PRIMARY 0.2c
gmt set MAP_LABEL_OFFSET 0.1c
gmt set MAP_TICK_LENGTH_PRIMARY -0.08c
gmt set MAP_FRAME_WIDTH 0.001c
gmt set MAP_FRAME_PEN thinner,black
gmt set MAP_GRID_PEN_PRIMARY 0.2p,darkgray,-
gmt set FONT_TITLE 20p,Helvetica,black
gmt set MAP_TITLE_OFFSET 0.1c

PS="${pngflname}.ps"
J=X18c/3c
Re=${stime}/${etime}/${min_e}/${max_e}
Rn=${stime}/${etime}/${min_n}/${max_n}
Ru=${stime}/${etime}/${min_u}/${max_u}
Rdop=${stime}/${etime}/${min_dop}/${max_dop}
Rnum=${stime}/${etime}/${min_sat}/${max_sat}

gmt psxy -J$J -R${Ru} -T -K > $PS

#E N U
gmt psbasemap -J$J -R${Re} -BWsen -Bxaf -Byaf+l"East (cm)" -Y11c -K -O >> $PS
cat enutmp | awk '{print $2,$3*100}' | gmt psxy -J$J -R${Renu} -Fas  -W0.03c -K -O >> $PS
echo East  ${rmse} cm | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PS

gmt psbasemap -J$J -R${Rn} -BWsen -Bxaf -Byaf+l"North (cm)" -Y-3.2c -K -O >> $PS
cat enutmp | awk '{print $2,$4*100}' | gmt psxy -J$J -R${Renu} -Fas  -W0.03c -K -O >> $PS
echo North  ${rmsn} cm | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PS

gmt psbasemap -J$J -R${Ru} -BWsen -Bxaf -Byaf+l"Up (cm)" -Y-3.2c -K -O >> $PS
cat enutmp | awk '{print $2,$5*100}' | gmt psxy -J$J -R${Renu} -Fas  -W0.03c -K -O >> $PS
echo Up  ${rmsu} cm | gmt pstext -R -J -O -K -F+cTR+jTR+f12p -Dj0.05i  >>$PS

#DOP Satnum

gmt psbasemap -J$J -R${Rdop} -BWSn -Bxaf+l"Seconds of day (s)" -Byaf+l"DOP" -Y-3.2c -K -O >> $PS
gmt psbasemap -J$J -R${Rnum} -BE -Byaf+l"Satellite number" -K -O >> $PS
cat enutmp | awk '{print $2,$7}' | gmt psxy -J$J -R${Rdop} -Fas -Wred -W0.03c -K -O >> $PS
cat enutmp | awk '{print $2,$6}' | gmt psxy -J$J -R${Rnum} -Fas -Wblue -W0.03c -K -O >> $PS
> legend.txt
echo "N 2" >> legend.txt
echo "S -3.3c c 0.1c red 1p,red,solid -3.1c PDOP" >> legend.txt
echo "S -2.5c c 0.1c blue 1p,blue,solid -2.3c Satellite number" >> legend.txt
cat legend.txt | gmt pslegend -J -R -DjTR+w1.5c+l0.8+o0.1c/0.1c  -K -O >> $PS

gmt psxy -J$J -R${Rnum} -T -O >> $PS

#gmt clear history
rm *.history
gmt psconvert -Tg -A -P -Z $PS

rm -f gmt.*  *.eps *.ps
rm -f enutmp legend.txt
