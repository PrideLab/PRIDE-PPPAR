#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Test PRIDE PPP-AR 2                                             ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 2.1                                                         ##
##                                                                           ##
##  DATE   : Sep-1, 2021                                                     ##
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
source ${HOME}/.bash_profile
lsq > /dev/null 2>&1
if [ $? -eq 127 ]; then  # command not found
    printf "${RED}error:${NC} PRIDE PPP-AR 2:lsq not found\n"
    printf "${RED}error:${NC} PRIDE PPP-AR 2 testing failed\n"; exit
fi

wk_dir=$(pwd)                  # working directory
# Generate control file
sed -i "s#/home/username/path-to-data#${wk_dir}#" config_template_daily > /dev/null 2>&1
sed -i "s#/home/username/path-to-table#${wk_dir}/..#" config_template_daily > /dev/null 2>&1
sed -i "s#/home/username/path-to-product#${wk_dir}#" config_template_daily > /dev/null 2>&1
sed -i "s#/home/username/path-to-data#${wk_dir}#" config_template_hourly > /dev/null 2>&1
sed -i "s#/home/username/path-to-table#${wk_dir}/..#" config_template_hourly > /dev/null 2>&1
sed -i "s#/home/username/path-to-product#${wk_dir}#" config_template_hourly > /dev/null 2>&1
sed -i "s#/home/username/path-to-data#${wk_dir}#" config_template_mobile > /dev/null 2>&1
sed -i "s#/home/username/path-to-table#${wk_dir}/..#" config_template_mobile > /dev/null 2>&1
sed -i "s#/home/username/path-to-product#${wk_dir}#" config_template_mobile > /dev/null 2>&1
sed -i "s#/home/username/path-to-data#${wk_dir}#" config_template_seismic > /dev/null 2>&1
sed -i "s#/home/username/path-to-table#${wk_dir}/..#" config_template_seismic > /dev/null 2>&1
sed -i "s#/home/username/path-to-product#${wk_dir}#" config_template_seismic > /dev/null 2>&1

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

echo -e "${BLUE}(1) static daily fixed"
pdp3 config_template_daily 20200101 20200101 abmf s 30 y
MvDir 2020/001 ./results/static-24h-fixed

echo -e "\n${BLUE}(2) kinematic daily fixed"
pdp3 config_template_daily 20210808 20210808 bako k 30 y
MvDir 2021/220 ./results/kinematic-24h-fixed

echo -e "\n${BLUE}(3) kinematic fixed (LAMBDA) 1h"
pdp3 config_template_hourly 20210729 20210729 ccj2 k 30 y
MvDir 2021/210 ./results/kinematic-1h-fixed-LAMBDA

echo -e "\n${BLUE}(4) high-rate fixed (LAMBDA) 1h"
pdp3 config_template_seismic 20210729 20210729 ac12 k 1 y
MvDir 2021/210 ./results/highrate-1h-fixed-LAMBDA

echo -e "\n${BLUE}(5) troposphere daily"
pdp3 config_template_daily 20200103 20200103 abpo f 30 y
MvDir 2020/003 ./results/tropo-24h-fixed

rm -rf 2020 2021

# Output
printf "${BLUE}::${NC} computation results are put in %s\n" ./results/
printf "${BLUE}::${NC} reference results are in %s\n" ./results_ref/
