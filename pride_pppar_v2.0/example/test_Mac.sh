#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Test PRIDE PPP-AR 2                                             ##
##                                                                           ##
##  AUTHOR : Yuanxin Pan    yxpan@whu.edu.cn                                 ##
##           For GPS-only PPP                                                ##
##                                                                           ##
##  MODIFIED BY: Songfeng Yang   sfyang@whu.edu.cn                           ##
##           For G/R/E/C/J PPP                                               ##
##                                                                           ##
##  VERSION: ver 2.0                                                         ##
##                                                                           ##
##  DATE   : Mar-24, 2021                                                    ##
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


RED='\033[0;31m'
BLUE='\033[1;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check installation
source ${HOME}/.bashrc
lsq > /dev/null 2>&1
if [ $? -eq 127 ]; then  # command not found
    printf "${RED}error:${NC} PRIDE PPP-AR 2:lsq not found\n"
    printf "${RED}error:${NC} PRIDE PPP-AR 2 testing failed\n"; exit
fi

wk_dir=$(pwd)                  # working directory
# Generate control file
config=$(basename `mktemp -u`)
config=${config/tmp/config}    # configuration file
echo "## Session configure" > $config
echo "Interval = 30"       >> $config
echo "Session time    = -YYYY- -MM- -DD- 00 00 00 86360" >> $config
echo "Rinex directory = ${wk_dir}/data/-YEAR-/-DOY-" >> $config
echo "Sp3 directory   = ${wk_dir}/products" >> $config
echo "Table directory = ${wk_dir}/../table" >> $config
cat config_partial >> $config

MvDir() {
    local src="$1"
    local dst="$2"
    local tmp=$(basename `mktemp -u`)
    tmp=${tmp/tmp/}
    [ -d "$dst" ] && mv "$dst" "${dst}${tmp}" && \
        echo -e "${BLUE}::${NC} mv existing $dst to ${dst}${tmp}"
    mv "$src" "$dst"
}

mkdir -p products results
# Computation
echo -e "(1) static float"
pride_pppar ${config} 20200101 20200101 N
MvDir 2020/001 ./results/static-float

echo -e "\n(2) static fixed"
pride_pppar ${config} 20200101 20200101 Y
MvDir 2020/001 ./results/static-fixed

sed -i '' '/^ abmf/ s/abmf S/abmf K/' ${config}
echo -e "\n(3) kinematic float"
pride_pppar ${config} 20200101 20200101 N
MvDir 2020/001 ./results/kinematic-float

echo -e "\n(4) kinematic fixed"
pride_pppar ${config} 20200101 20200101 Y
MvDir 2020/001 ./results/kinematic-fixed

sed -i '' '/Session time/s/86360/3600/' ${config}
echo -e "\n(5) kinematic float 1h"
pride_pppar ${config} 20200101 20200101 N
MvDir 2020/001 ./results/kinematic-1h-float

sed -i '' '/Remove bias/s/YES/NO/' ${config}
sed -i '' '/Ambiguity fixing/s/ROUNDING/LAMBDA/' ${config}
echo -e "\n(6) kinematic fixed (LAMBDA) 1h"
pride_pppar ${config} 20200101 20200101 Y
MvDir 2020/001 ./results/kinematic-1h-fixed-LAMBDA

sed -i '' '/Session time/s/3600/86360/' ${config}
sed -i '' '/Interval/s/30/0.5/' ${config}
sed -i '' '/ZTD model/s/PWC:60/STO/' ${config}
sed -i '' '/HTG model/s/PWC:720/STO/' ${config}
sed -i '' '/^ abmf/ s/abmf K/rov1 K/' ${config}
echo -e "\n(7) high-dynamic"
pride_pppar ${config} 20081005 20081005 N
MvDir 2008/279 ./results/high-dynamic

rm -rf 2020 2008

# Output
printf "${BLUE}::${NC} computation results are put in %s\n" ./results/
printf "${BLUE}::${NC} reference results are in %s\n" ./results_ref/
