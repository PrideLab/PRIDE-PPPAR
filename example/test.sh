#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Test PRIDE PPP-AR                                               ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 3.0                                                         ##
##                                                                           ##
##  DATE   : Sept-13, 2023                                                   ##
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


RED='\033[0;31m'
BLUE='\033[1;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check installation
lsq > /dev/null 2>&1
if [ $? -eq 127 ]; then  # command not found
    printf "${RED}error:${NC} PRIDE PPP-AR 3 lsq not found\n"
    printf "${RED}error:${NC} PRIDE PPP-AR 3 testing failed\n"; exit
fi

MvDir() {
    local src="$1"
    local dst="$2"
    local tmp=$(basename `mktemp -u`)
    tmp=${tmp/tmp/}
    [ -d "$dst" ] && mv "$dst" "${dst}${tmp}" && \
        echo -e "${BLUE}::${NC} mv existing $dst to ${dst}${tmp}"
    mv "$src" "$dst"
}

mkdir -p results

# Computation
echo -e "${BLUE}(1) static daily fixed"
pdp3 -m S ./data/2020/001/abmf0010.20o
MvDir 2020/001 ./results/static-24h-fixed

echo -e "\n${BLUE}(2) kinematic daily fixed"
pdp3 ./data/2021/220/BAKO00IDN_R_20212200000_01D_30S_MO.rnx
MvDir 2021/220 ./results/kinematic-24h-fixed

echo -e "\n${BLUE}(3) kinematic fixed (LAMBDA) 1h"
pdp3 -z P ./data/2021/210/ccj22100.21o
MvDir 2021/210 ./results/kinematic-1h-fixed-LAMBDA

echo -e "\n${BLUE}(4) high-rate fixed (LAMBDA) 1h"
pdp3 ./data/2021/210/ac122100.21o
MvDir 2021/210 ./results/highrate-1h-fixed-LAMBDA

echo -e "\n${BLUE}(5) troposphere daily"
pdp3 -m F ./data/2020/003/abpo0030.20o
MvDir 2020/003 ./results/tropo-24h-fixed

rm -rf 2020 2021

# Output
printf "${BLUE}::${NC} computation results are put in %s\n" ./results/
printf "${BLUE}::${NC} reference results are in %s\n" ./results_ref/
