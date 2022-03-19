#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Convert weekly combination of IGS daily combined solutions      ##
##           to the required format of PRIDE PPP-AR (sit.xyz)                ##
##                                                                           ##
##  AUTHOR : PRIDE LAB      pride@whu.edu.cn                                 ##
##                                                                           ##
##  VERSION: ver 2.1                                                         ##
##                                                                           ##
##  DATE   : Sept-24, 2021                                                   ##
##                                                                           ##
##              @ GNSS RESEARCH CENTER, WUHAN UNIVERSITY, 2021               ##
##                                                                           ##
##    Copyright (C) 2021 by Wuhan University                                 ##
##                                                                           ##
##    This program is free software: you can redistribute it and/or modify   ##
##    it under the terms of the GNU General Public License (version 3) as    ##
##    published by the Free Software Foundation.                             ##
##                                                                           ##
##    This program is distributed in the hope that it will be useful,        ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of         ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          ##
##    GNU General Public License (version 3) for more details.               ##
##                                                                           ##
##    You should have received a copy of the GNU General Public License      ##
##    along with this program.  If not, see <https://www.gnu.org/licenses/>. ##
##                                                                           ##
###############################################################################

main()
{

    CheckCmdArgs "$@" || exit 1
	
    local site=`echo $1 | awk '{print toupper($1)}'`
    local year=$2
    local yr=`echo $2 | cut -c3-4`
    local mm=$3
    local dd=$4
    local dow=`ymd2wkdow $year $mm $dd | awk '{printf("%d\n",$2)}'`
    local week=`ymd2wkdow $year $mm $dd | awk '{printf("%d\n",$1)}'`
    local args="-nv -t 3 --connect-timeout=10 --read-timeout=60"
    url=ftp://gssc.esa.int/gnss

	if [ ! -e igs${yr}P${week}.ssc ]; then
      wget ${args} ${url}/products/${week}/igs${yr}P${week}.ssc.Z \
        && gunzip igs${yr}P${week}.ssc.Z
	fi
    [ ! -e igs${yr}P${week}.ssc ] && echo "ERROR: No igs${yr}P${week}.ssc" && return 1
    
    awk -v sit=$site 'BEGIN{fg=0;x=0.0;y=0.0;z=0.0;sigx=0.0;sigy=0.0;sigz=0.0;snam=" "}\
              {\
                if($1=="+SOLUTION/ESTIMATE"){fg=1};if($1=="-SOLUTION/ESTIMATE"){fg=0};\
                if(fg==1)\
                {\
                  if($2=="STAX"){snam=$3;x=$9;sigx=$10};\
                  if($2=="STAY"){y=$9;sigy=$10};\
                  if($2=="STAZ"&&$3==sit)\
                  {\
                    z=$9;sigz=$10;printf(" %25.6f %25.6f %25.6f %25.6f %25.6f %25.6f\n",x,y,z,sigx,sigy,sigz);\
                    snam=" ";x=0.0;y=0.0;z=0.0;sigx=0.0;sigy=0.0;sigz=0.0;\
                  };\
               }\
             }' igs${yr}P${week}.ssc

}

CheckCmdArgs() { # purpose: chech whether command line arguments are right
                 # usage  : CheckCmdArgs "$@"
    if [ $# -ne 4 ]; then
        echo "Usage : snxsit.sh SITENAME YYYY MM DD"
        return 1
    else
	local ymd=$2$3$4
    local year=$2
	local mm=$3
    local dd=$4

    if [ ${#ymd} -ne 8 ]; then
       echo -e "ERROR: Bad date format: $ymd"
       echo -e "Date format: YYYY MM DD"
       return 1
	fi
	echo "$ymd" | [ -n "`sed -n '/^[0-9][0-9]*$/p'`" ] 
	if [ $? -ne 0 ]; then
	  echo -e "ERROR: Input date is not a number: $ymd"
	  return 1
	elif [ $year -le 1980 ]; then
	  echo -e "ERROR: Input year is wrong $year"
	  return 1
	elif [ $mm -le 0 -o $mm -gt 12 ]; then
	  echo -e "ERROR: Input month is wrong : $mm"
	  return 1
	elif [ $dd -le 0 -o $dd -gt 31 ]; then
	  echo -e "ERROR: Input day is wrong : $dd"
	  return 1
	fi
    fi
}

ymd2mjd()
{
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
    #local mjd=$(bc <<< "$year*365.25 - $year*365.25 % 1 - 679006")
    #mjd=$(bc <<< "($mjd + (30.6001*($mon+1))/1 + 2 - $year/100 + $year/400 + $day)/1")
    echo $mjd
}
ymd2wkdow()
{
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
######################################################################
##                               Entry                              ##
######################################################################
main "$@"
