#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Prepare LEO RINEX observation files and products                ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 3.1.0                                                        ##
##                                                                           ##
##  DATE   : Jan-04, 2025                                                    ##
##                                                                           ##
##              @ GNSS RESEARCH CENTER, WUHAN UNIVERSITY, 2025               ##
##                                                                           ##
##    Copyright (C) 2025 by Wuhan University                                 ##
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


main(){

if [ $# -ne 2 -a $# -ne 1 ]; then
    echo "#usage  : ./plotmhm.sh mhm_filename png_filename"
    echo "#example: ./plotmhm.sh mhm_rov1 png_mhm_rov1"
    echo "#if no png_filename, default mhm_filename_P or mhm_filename_L"
    exit
fi

mhmflname=$1
if [ $# -eq 2 ]; then
    pngflname=$2
elif [ $# -eq 1 ]; then
    pngflname=$1
#    pngflname=`echo $* | awk '{print $1}'`
fi

site=""
daystart=""
dayend=""
sysnum=0
first_daylist=false
inside_part=false
while IFS= read -r line; do
    # read header
    if [[ "$line" == *"STATION"* ]]; then
        site="${line:0:4}"
        site=$(echo $site | tr '[:lower:]' '[:upper:]')
    fi
    if [[ "$line" == *"MODELING DURATION (DAY)"* ]]; then
        daynum=`echo $line | awk '{print $1}'`
    fi
    if [[ "$line" == *"DAY LIST"* ]]; then
        if [ "$first_daylist" = false ]; then
            daystart=`echo $line | awk '{print $1}'`
            first_daylist=true
        fi
        dayend=`echo $line | awk '{print $(NF-2)}'`
    fi
    
    # prepare data
    if [[ "$line" == *"START OF"* ]]; then
        inside_part=true
        comnum=0
        tmpdata="tmp_datamhm"
        > $tmpdata    
        sysarray=($(echo "$line" | sed 's/START OF//g' | tr -s ' ' '\n'))
        for i in "${!sysarray[@]}"; do
            if [[ $((i % 2)) -eq 0 && "${sysarray[$i]}" != "QZS" ]]; then
                sys="${sysarray[$i]}"
                sysname+=("$sys")
                ((sysnum=sysnum+1))
                ((comnum=comnum+1))
                tmpdataP="tmp_${sys}P"
                tmpdataL="tmp_${sys}L"
                > $tmpdataP
                > $tmpdataL
            fi
        done
    elif [[ "$line" == *"END OF"* ]]; then
        if $inside_part; then
            ((j=$sysnum-$comnum))
            for ((i=$j;i<$sysnum;i++)); do                
                sys="${sysname[$i]}"
                awk 'NR > 1 {print $1, 90 - $2, $3}' "${tmpdata}" > "tmp_${sys}P"
                awk 'NR > 1 {print $1, 90 - $2, $7}' "${tmpdata}" > "tmp_${sys}L"
            done
            inside_part=false
        fi
    elif $inside_part; then
        echo $line >> ${tmpdata}
    fi        
done < "$mhmflname"

# title
for i in "${!sysname[@]}"; do
    sys=${sysname[$i]}
    if [ "$sys" == "GPS" ]; then
        axtitle[$i]="GPS"
    elif [ "$sys" == "GAL" ]; then
        axtitle[$i]="Galileo"
    elif [ "$sys" == "BDS" ]; then
        axtitle[$i]="BDS"
    elif [ "$sys" == "GLO" ]; then
        axtitle[$i]="Glonass"
    else
        axtitle[$i]=$sys
    fi
done

mode[0]="P"
mode[1]="L"
cptrange[0]="-1.5/1.5"
cptrange[1]="-0.03/0.03"
scale[0]="0.5"
scale[1]="0.015"
title[0]="Pseudorange Multipath Correction (meter) "
title[1]="Carrier Phase Multipath Correction (meter) "

if [ $daynum -eq 1 ]; then
    daylist="${site} ${daystart}"
elif [ $daynum -gt 1 ]; then
    daylist="${site} ${daystart}-${dayend}"
fi

# sysnum = 1 , one figure
if [ $sysnum -eq 1 ]; then
for ((i=0;i<2;i++)); do
# parameters
gmt set FORMAT_GEO_MAP +D 
gmt set FONT Helvetica
gmt set FONT_TITLE 7p,0,black
gmt set FONT_LABEL 8p,0,black
gmt set FONT_ANNOT 4p,0,black
gmt set MAP_TITLE_OFFSET 12p
gmt set MAP_GRID_PEN default,grey
gmt set MAP_ANNOT_ORTHO we
gmt set MAP_FRAME_PEN 0.3p,black
gmt set MAP_TICK_LENGTH 0p/0p

PS="${pngflname}_${mode[$i]}.ps"

# Azimuth axes
gmt pstext -R0/10/0/4 -JX10c/4c -F+f4p,0,black -P -K > $PS << EOF
1.85 3.25 0\232
3.12 2.57 60\232
3.12 1 120\232
1.85 0.35 180\232
0.5 1 240\232
0.5 2.57 300\232
EOF

# plotmhm
gmt makecpt -Cjet -T${cptrange[$i]} > mhmjet.cpt
gmt psxy tmp_${sysname[0]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X0.5c -Y0.5c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[0]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black+jMC -K -O >> $PS << EOF   
0 30 60\232
0 60 30\232
EOF

# cocorbar
gmt psscale -R0/10/0/4 -JX10c/4c -DjTC+w5.5c/0.15c+o-3.7c/4.5c -Bxa${scale[$i]}f${scale[$i]}+l"${title[$i]} ${daylist}" -Cmhmjet.cpt --MAP_TICK_LENGTH=1.5p/1p -O >> $PS   


# PNG
gmt psconvert -A+S0.7 -E700 -Tg $PS

rm -f *.history 
rm -f gmt.*  *.eps *.ps
rm -f mhmjet.cpt
rm -f $tmpdata
rm -f tmp_${sysname[0]}${mode[$i]}
done



# sysnum = 2 , two figures
elif [ $sysnum -eq 2 ]; then
for ((i=0;i<2;i++)); do
# parameters
gmt set FORMAT_GEO_MAP +D 
gmt set FONT Helvetica
gmt set FONT_TITLE 7p,0,black
gmt set FONT_LABEL 8p,0,black
gmt set FONT_ANNOT 4p,0,black
gmt set MAP_TITLE_OFFSET 12p
gmt set MAP_GRID_PEN default,grey  
gmt set MAP_ANNOT_ORTHO we
gmt set MAP_FRAME_PEN 0.3p,black
gmt set MAP_TICK_LENGTH 0p/0p

PS="${pngflname}_${mode[$i]}.ps"

# Azimuth axes
gmt pstext -R0/10/0/4 -JX10c/4c -F+f4p,0,black -P -K > $PS << EOF
1.85 3.25 0\232
3.12 2.57 60\232
3.12 1 120\232
1.85 0.35 180\232
0.5 1 240\232
0.5 2.57 300\232
5.05 3.25 0\232
6.32 2.57 60\232
6.32 1 120\232
5.05 0.35 180\232
3.7 1 240\232
3.7 2.57 300\232
EOF

# plotmhm
gmt makecpt -Cjet -T${cptrange[$i]} > mhmjet.cpt
gmt psxy tmp_${sysname[0]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X0.5c -Y0.5c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[0]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF 
0 30 60\232
0 60 30\232
EOF

gmt psxy tmp_${sysname[1]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X3.2c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[1]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF
0 30 60\232
0 60 30\232
EOF

# cocorbar
gmt psscale -R0/10/0/4 -JX10c/4c -DjTC+w6c/0.15c+o-5.3c/4.5c -Bxa${scale[$i]}f${scale[$i]}+l"${title[$i]} ${daylist}" -Cmhmjet.cpt --MAP_TICK_LENGTH=1.5p/1p -O >> $PS     

# PNG
gmt psconvert -A -E650 -Tg $PS


rm -f gmt.*
rm -f *.ps
rm -f mhmjet.cpt
rm -f $tmpdata
rm -f tmp_${sysname[0]}${mode[$i]}
rm -f tmp_${sysname[1]}${mode[$i]}
done



# sysnum = 3 , three figures
elif [ $sysnum -eq 3 ]; then
for ((i=0;i<2;i++)); do

gmt set FORMAT_GEO_MAP +D 
gmt set FONT Helvetica
gmt set FONT_TITLE 7p,0,black
gmt set FONT_LABEL 8p,0,black
gmt set FONT_ANNOT 4p,0,black
gmt set MAP_TITLE_OFFSET 12p
gmt set MAP_GRID_PEN default,grey
gmt set MAP_ANNOT_ORTHO we
gmt set MAP_FRAME_PEN 0.3p,black
gmt set MAP_TICK_LENGTH 0p/0p

PS="${pngflname}_${mode[$i]}.ps"

# Azimuth axes
gmt pstext -R0/10/0/4 -JX10c/4c -F+f4p,0,black -P -K > $PS << EOF
1.85 3.25 0\232
3.12 2.57 60\232
3.12 1 120\232
1.85 0.35 180\232
0.5 1 240\232
0.5 2.57 300\232
5.05 3.25 0\232
6.32 2.57 60\232
6.32 1 120\232
5.05 0.35 180\232
3.7 1 240\232
3.7 2.57 300\232
8.25 3.25 0\232
9.52 2.57 60\232
9.52 1 120\232
8.25 0.35 180\232
6.9 1 240\232
6.9 2.57 300\232
EOF

# plotmhm
gmt makecpt -Cjet -T${cptrange[$i]} > mhmjet.cpt
gmt psxy tmp_${sysname[0]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X0.5c -Y0.5c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[0]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF   
0 30 60\232
0 60 30\232
EOF

gmt psxy tmp_${sysname[1]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X3.2c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[1]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF
0 30 60\232
0 60 30\232
EOF

gmt psxy tmp_${sysname[2]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X3.2c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[2]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF
0 30 60\232
0 60 30\232
EOF

# colorbar
gmt psscale -R0/10/0/4 -JX10c/4c -DjTC+w7c/0.15c+o-6.9c/4.5c -Bxa${scale[$i]}f${scale[$i]}+l"${title[$i]} ${daylist}" -Cmhmjet.cpt --MAP_TICK_LENGTH=1.5p/1p -O >> $PS     

# PNG
gmt psconvert -A -E500 -Tg $PS

rm -f gmt.*
rm -f *.ps
rm -f mhmjet.cpt
rm -f $tmpdata
rm -f tmp_${sysname[0]}${mode[$i]}
rm -f tmp_${sysname[1]}${mode[$i]}
rm -f tmp_${sysname[2]}${mode[$i]}
done


# sysnum = 4 , four figures
elif [ $sysnum -eq 4 ]; then
for ((i=0;i<2;i++)); do

gmt set FORMAT_GEO_MAP +D 
gmt set FONT Helvetica
gmt set FONT_TITLE 7p,0,black
gmt set FONT_LABEL 8p,0,black
gmt set FONT_ANNOT 4p,0,black
gmt set MAP_TITLE_OFFSET 12p
gmt set MAP_GRID_PEN default,grey 
gmt set MAP_ANNOT_ORTHO we
gmt set MAP_FRAME_PEN 0.3p,black
gmt set MAP_TICK_LENGTH 0p/0p

PS="${pngflname}_${mode[$i]}.ps"

# Azimuth axes
gmt pstext -R0/14/0/4 -JX14c/4c -F+f4p,0,black -P -K > $PS << EOF
1.85 3.25 0\232
3.12 2.57 60\232
3.12 1 120\232
1.85 0.35 180\232
0.5 1 240\232
0.5 2.57 300\232
5.05 3.25 0\232
6.32 2.57 60\232
6.32 1 120\232
5.05 0.35 180\232
3.7 1 240\232
3.7 2.57 300\232
8.25 3.25 0\232
9.52 2.57 60\232
9.52 1 120\232
8.25 0.35 180\232
6.9 1 240\232
6.9 2.57 300\232
11.45 3.25 0\232
12.72 2.57 60\232
12.72 1 120\232
11.45 0.35 180\232
10.1 1 240\232
10.1 2.57 300\232
EOF

# plotmhm
gmt makecpt -Cjet -T${cptrange[$i]} > mhmjet.cpt
gmt psxy tmp_${sysname[0]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X0.5c -Y0.5c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[0]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF      # -D文本偏移量
0 30 60\232
0 60 30\232
EOF

gmt psxy tmp_${sysname[1]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X3.2c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[1]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF
0 30 60\232
0 60 30\232
EOF

gmt psxy tmp_${sysname[2]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X3.2c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[2]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF
0 30 60\232
0 60 30\232
EOF

gmt psxy tmp_${sysname[3]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X3.2c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[3]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF
0 30 60\232
0 60 30\232
EOF

# colorbar
gmt psscale -R0/14/0/4 -JX14c/4c -DjTC+w7c/0.15c+o-10.5c/4.5c -Bxa${scale[$i]}f${scale[$i]}+l"${title[$i]} ${daylist}" -Cmhmjet.cpt --MAP_TICK_LENGTH=1.5p/1p -O >> $PS     

# PNG
gmt psconvert -A -E500 -Tg $PS


rm -rf gmt.*
rm -rf *.ps
rm -rf mhmjet.cpt
rm -f $tmpdata
rm -f tmp_${sysname[0]}${mode[$i]}
rm -f tmp_${sysname[1]}${mode[$i]}
rm -f tmp_${sysname[2]}${mode[$i]}
rm -f tmp_${sysname[3]}${mode[$i]}
done


# sysnum > 4 , more than 4 figures
elif [ $sysnum -gt 4 ]; then
for ((i=0;i<2;i++)); do

gmt set FORMAT_GEO_MAP +D 
gmt set FONT Helvetica
gmt set FONT_TITLE 7p,0,black
gmt set FONT_LABEL 8p,0,black
gmt set FONT_ANNOT 4p,0,black
gmt set MAP_TITLE_OFFSET 12p
gmt set MAP_GRID_PEN default,grey
gmt set MAP_ANNOT_ORTHO we
gmt set MAP_FRAME_PEN 0.3p,black
gmt set MAP_TICK_LENGTH 0p/0p

PS="${pngflname}_${mode[$i]}.ps"
pngwidth=$(echo "10 + 4 * ($sysnum - 3)" | bc)
R="0/${pngwidth}/0/4"
J="X${pngwidth}c/4c"


# Azimuth axes
tmpcoor="tmp_coordata"
> $tmpcoor
coordata=(
"1.85 3.25 0\232"
"3.12 2.57 60\232"
"3.12 1 120\232"
"1.85 0.35 180\232"
"0.5 1 240\232"
"0.5 2.57 300\232")

for ((j=0;j<$sysnum;j++)); do
    for line in "${coordata[@]}"; do
        x_coor=$(echo "$line" | cut -d' ' -f1)
        rest_coor=$(echo "$line" | cut -d' ' -f2-)
        new_coor=$(echo "$x_coor + 3.2 * $j" | bc)
        echo "$new_coor $rest_coor" >> $tmpcoor
    done
done
gmt pstext $tmpcoor -R$R -J$J -F+f4p,0,black -P -K > $PS


# plotmhm
gmt makecpt -Cjet -T${cptrange[$i]} > mhmjet.cpt
gmt psxy tmp_${sysname[0]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X0.5c -Y0.5c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[0]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF   
0 30 60\232
0 60 30\232
EOF

for ((j=1;j<$sysnum;j++)); do
gmt psxy tmp_${sysname[$j]}${mode[$i]} -JPa2.6c -R0/360/0/90 -Sc0.06c -Cmhmjet.cpt -X3.2c -K -O >> $PS
gmt psbasemap -JPa2.6c -R0/360/0/90 -Bxf60g60 -Byf30g30 -B+t"${axtitle[$j]}" -K -O >> $PS
gmt pstext -JPa2.6c -R0/360/0/90 -D0.17c/0.1c -F+f4p,0,black -K -O >> $PS << EOF
0 30 60\232
0 60 30\232
EOF
done


# colorbar
barmove=$(echo "6.9 + 3.6 * ($sysnum - 3)" | bc)
gmt psscale -R$R -J$J -DjTC+w7c/0.15c+o-${barmove}c/4.5c -Bxa${scale[$i]}f${scale[$i]}+l"${title[$i]} ${daylist}" -Cmhmjet.cpt --MAP_TICK_LENGTH=1.5p/1p -O >> $PS     

# PNG
gmt psconvert -A -E500 -Tg $PS



rm -f gmt.*
rm -f *.ps
rm -f mhmjet.cpt
rm -f $tmpdata
rm -f $tmpcoor
for ((j=0;j<$sysnum;j++)); do
rm -f tmp_${sysname[$j]}${mode[$i]}
done

done
if [ $sysnum -gt 5 ]; then
echo "Warning: Too many subplots are drawn and may exceed the maximum width of the paper"
fi

fi































}
main "$@"
