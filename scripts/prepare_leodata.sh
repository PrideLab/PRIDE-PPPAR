#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Prepare LEO RINEX observation files and products                ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 3.0                                                         ##
##                                                                           ##
##  DATE   : Sept-12, 2023                                                   ##
##                                                                           ##
##              @ GNSS RESEARCH CENTER, WUHAN UNIVERSITY, 2023               ##
##                                                                           ##
##    Copyright (C) 2023 by Wuhan University                                 ##
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

readonly NC='\033[0m'
readonly RED='\033[0;31m'
readonly BLUE='\033[1;34m'
readonly CYAN='\033[0;36m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'

readonly MSGERR="${RED}error:$NC"
readonly MSGWAR="${YELLOW}warning:$NC"
readonly MSGINF="${BLUE}::$NC"
readonly MSGSTA="${BLUE}===>$NC"

readonly LINSTALL=NO
readonly softpath="${HOME}/software"

pathfile() {
   local pfile
   if [ "$OS" == "Darwin" ]; then
      echo ${HOME}/.bash_profile
   else
      echo ${HOME}/.bashrc
   fi
}

main() {
   if [ $# -ne 2 ]; then
      >&2 echo "Usage      : prepare_leodata.sh leo_code date"
      >&2 echo "  leo_name : LEO mission name      select one from: grace/grace-fo"
      >&2 echo "  date     : date                  input as year/month/day or year/doy(day of year)"
      >&2 echo "Example    : $0 grace-fo 2023/1"
      exit 1
   fi
   
   local leo_code=$(echo $1 | tr 'a-z' 'A-Z')
   local date="$2"
   
   [ ${LINSTALL} == "YES" ] && CheckExecutables
  
   local time=($(echo $date | tr '/:-' ' '))
   case ${#time[@]} in
   2 ) ymd=($(ydoy2ymd ${time[@]}  | awk '{printf("%04d %02d %02d\n",$1,$2,$3)}'));;
   3 ) ymd=($(echo    "${time[@]}" | awk '{printf("%04d %02d %02d\n",$1,$2,$3)}')) ;;
   * ) >&2 echo -e "$MSGERR downleo.sh: invalid time format: $date" && exit 1
   esac
   mjd=$(ymd2mjd ${ymd[@]})
   
   [ ! -d data/${ymd[0]} ] && mkdir -p data/${ymd[0]}
   [ ! -d ${ymd[0]}/product/leo ] && mkdir -p ${ymd[0]}/product/leo
   otherfiledir="${leo_code}-files"
   [ ! -d ${otherfiledir}/${ymd[0]} ] && mkdir -p ${otherfiledir}/${ymd[0]}

   case ${leo_code} in
   GRACE ) PrepareGRACE $mjd ;;
   GRACE-FO ) PrepareGRACEFO $mjd ;;
   * ) >&2 echo -e "$MSGERR invalid LEO id: ${leo_code}" && exit 1
   esac
}

CheckExecutables() {
   ## GRACE: Bin2AsciiLevel1.e or gps1x2rnx.e
   for exceu in "Bin2AsciiLevel1.e" "gps1x2rnx.e"; do
       if ! which $exceu > /dev/null 2>&1; then
           echo -e -n "$MSGWAR "
           read -p "File format transform tool $exceu for GRACE not found, install or not [Y/N]:" lgrace
           if [ ${#lgrace} -ge 1 ] && ( [ ${lgrace:0:1} == "y" ] || [ ${lgrace:0:1} == "Y" ] ); then
               InstallSoftGRACE
           fi
       fi
   done
}

InstallSoftGRACE() {
   local pdp3dir="${HOME}/.PRIDE_PPPAR_BIN/"
   local cwd=$(pwd)
   local GraceReadSW="GraceReadSW_L1_2010-03-31.tar.gz"
   local GraceReadurl="ftp://isdcftp.gfz-potsdam.de/grace/SOFTWARE/${GraceReadSW}"
   local GraceReadPath="${softpath}/RELEASE_2010-03-31"
   local pathf=$(pathfile)
   
   [ ! -d "$softpath" ] && mkdir "$softpath"
   cd "$softpath"
   [ ! -f "$GraceReadSW" ] && WgetDownload "$GraceReadurl"
   
   tar -zxvf "$GraceReadSW"
   cd RELEASE_2010-03-31/
   sed -i "/GRACESystemPath/s# = .*# = \"${GraceReadPath}\";#" GRACEsyspath.h
   make 
   grep "^export PATH=${GraceReadPath}:\$PATH" $pathf > /dev/null 2>&1
   [ $? -ne 0 ] && echo "export PATH=${GraceReadPath}:\$PATH" >> $pathf
   source $pathf
   cd "$cwd"
}

WgetDownload() { # purpose : download a file with wget
                 # usage   : WgetDownload url
    local url="$1"
    local arg="-q -nv -nc -c -t 3 --connect-timeout=10 --read-timeout=60"
    [ -n "$url" ] || return 1

    wget --help | grep -q "\-\-show\-progress" && arg="$arg --show-progress"
    local cmd="wget $arg $url"
    echo "$cmd" | bash

    [ -e $(basename "$url") ] && return 0  || return 1
}

PrepareGRACEFO() {
    local mjd=$1
    local ydoy=($(mjd2ydoy $mjd))
    local ymd=($(ydoy2ymd ${ydoy[*]}))
    local yy=${ymd[0]:2:2}

    if [ ${mjd} -lt 58260 ]; then
        >&2 echo -e "$MSGERR Input date out of range: ${date}, should be after 2018/04/04"
        exit 1
    fi

    local grac="gracefo_1B_${ymd[0]}-${ymd[1]}-${ymd[2]}_RL04.ascii.noLRI"
    local grac_cmb="${grac}.tgz"
    local url="ftp://isdcftp.gfz-potsdam.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/${ymd[0]}/${grac_cmb}"
    
    if [ ! -e $grac_cmb ]; then
        WgetDownload "$url"
        if [ $? -ne 0 ]; then
            >&2 echo -e "$MSGERR GRACE data download failed: ${mjd}"
            exit 1
        fi
    fi
    local file_prefix=("GPS1B" "GNV1B" "SCA1B")
    for prefix in "${file_prefix[@]}"; do
        suffix="txt"
        [ $prefix == "GPS1B" ] && suffix="rnx"
        local file1="${prefix}_${ymd[0]}-${ymd[1]}-${ymd[2]}_C_04.${suffix}"
        local file2="${prefix}_${ymd[0]}-${ymd[1]}-${ymd[2]}_D_04.${suffix}"
        #echo $file $file2
        tar zxvf "$grac_cmb" $file1
        tar zxvf "$grac_cmb" $file2
        case $prefix in
        "GPS1B" )
            local rnx1="grac${ydoy[1]}0.${yy}o"
            local rnx2="grad${ydoy[1]}0.${yy}o"
            mv $file1 data/${ymd[0]}/$rnx1
            mv $file2 data/${ymd[0]}/$rnx2
            ;;
        "SCA1B" )
            local lat1="lat_${ydoy[0]}${ydoy[1]}_grac"
            local lat2="lat_${ydoy[0]}${ydoy[1]}_grad"
            lat2obx.py $file1 grac && mv $lat1 ${ymd[0]}/product/leo/
            lat2obx.py $file2 grad && mv $lat2 ${ymd[0]}/product/leo/
            rm -f $file1 $file2
            ;;
        "GNV1B" )
            local pso1="pso_${ydoy[0]}${ydoy[1]}_grac"
            local pso2="pso_${ydoy[0]}${ydoy[1]}_grad"
            mv $file1 ${ymd[0]}/product/leo/$pso1
            mv $file2 ${ymd[0]}/product/leo/$pso2
            ;;
        esac
    done
    mv "$grac_cmb" ${otherfiledir}/${ymd[0]}
}

PrepareGRACE() {
    local mjd=$1
    local ydoy=($(mjd2ydoy $mjd))
    local ymd=($(ydoy2ymd ${ydoy[*]}))
    local yy=${ymd[0]:2:2}
  
    if [ ${mjd} -lt 52368 -o ${mjd} -gt 58057 ]; then
        >&2 echo -e "$MSGERR Input date out of range: ${date}, should be in 2002/04/04-2017/10/31"
        exit 1
    fi

    local grac="grace_1B_${ymd[0]}-${ymd[1]}-${ymd[2]}_02"
    local grac_cmb="${grac}.tar.gz"
    local url="ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/JPL/INSTRUMENT/RL02/${ymd[0]}/${grac_cmb}"

    if [ ! -e $grac_cmb ]; then
        WgetDownload "$url"
        if [ $? -ne 0 ]; then
            >&2 echo -e "$MSGERR GRACE-FO data download failed: ${mjd}"
            exit 1
        fi
    fi
    local file_prefix=("GPS1B" "GNV1B" "SCA1B")
    for prefix in "${file_prefix[@]}"; do
        local file1="${prefix}_${ymd[0]}-${ymd[1]}-${ymd[2]}_A_02.dat"
        local file2="${prefix}_${ymd[0]}-${ymd[1]}-${ymd[2]}_B_02.dat"
        #echo $file $file2
        tar zxvf "$grac_cmb" $file1
        tar zxvf "$grac_cmb" $file2
        case $prefix in
        "GPS1B" )
            local rnx1="graa${ydoy[1]}0.${yy}o"
            local rnx2="grab${ydoy[1]}0.${yy}o"
            which gps1x2rnx.e > /dev/null 2>&1
            if [ $? -eq 0 ]; then
               gps1x2rnx.e  -gps1x $file1 -rnx $rnx1 && mv $rnx1 data/${ymd[0]}
               gps1x2rnx.e  -gps1x $file2 -rnx $rnx2 && mv $rnx2 data/${ymd[0]}
            fi 
            mv $file1 data/${ymd[0]}
            mv $file2 data/${ymd[0]}
            ;; 
        "SCA1B" )
            local lat1="lat_${ydoy[0]}${ydoy[1]}_graa"
            local lat2="lat_${ydoy[0]}${ydoy[1]}_grab"
            which Bin2AsciiLevel1.e > /dev/null 2>&1
            if [ $? -eq 0 ]; then
                Bin2AsciiLevel1.e -binfile $file1 -ascfile $lat1 && lat2obx.py $lat1 graa && mv $lat1 ${ymd[0]}/product/leo/
                Bin2AsciiLevel1.e -binfile $file2 -ascfile $lat2 && lat2obx.py $lat2 grab && mv $lat2 ${ymd[0]}/product/leo/
            fi
            mv $file1 ${ymd[0]}/product/leo
            mv $file2 ${ymd[0]}/product/leo
            ;; 
        "GNV1B" )
            local pso1="pso_${ydoy[0]}${ydoy[1]}_graa"
            local pso2="pso_${ydoy[0]}${ydoy[1]}_grab"
            which Bin2AsciiLevel1.e > /dev/null 2>&1
            if [ $? -eq 0 ]; then
                Bin2AsciiLevel1.e -binfile $file1 -ascfile $pso1 && mv $pso1 ${ymd[0]}/product/leo/
                Bin2AsciiLevel1.e -binfile $file2 -ascfile $pso2 && mv $pso2 ${ymd[0]}/product/leo/
            fi  
            mv $file1 ${ymd[0]}/product/leo
            mv $file2 ${ymd[0]}/product/leo
            ;; 
        esac
    done
    mv "$grac_cmb" ${otherfiledir}/${ymd[0]}
} 

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

wkdow2ydoy() {
    local week=$((10#$1))
    local dow=$2
    local mjd0=44243
    local mjd=$(($mjd0+1+$week*7+$dow))
    mjd2ydoy $mjd
}
ymd2wkdow() {
    local year=$1
    local mon=$2
    local day=$3
    local mjd0=44243
    local mjd=$(ymd2mjd $year $mon $day)
    local difmjd=$(($mjd-$mjd0-1))
    local week=$(($difmjd/7))
    local dow=$(($difmjd%7))
    printf "%04d %d\n" $week $dow
}

mjd2wkdow() {
    local mjd=$1
    local mjd0=44243
    local difmjd=$(($mjd-$mjd0-1))
    local week=$(($difmjd/7))
    local dow=$(($difmjd%7))
    printf "%04d %d\n" $week $dow
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
