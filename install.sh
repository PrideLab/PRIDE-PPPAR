#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Install PRIDE-PPPAR                                             ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 3.0                                                         ##
##                                                                           ##
##  DATE   : Feb-04, 2024                                                    ##
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

SYS=`uname`
RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check compiler
gfortran --version > /dev/null 2>&1
if [ $? -ne 0 ]; then
    printf "${RED}error:${NC} no compiler: gfortran\n"
    printf "${RED}error:${NC} PRIDE-PPPAR installation failed\n"; exit
fi
make --version > /dev/null 2>&1
if [ $? -ne 0 ]; then
    printf "${RED}error:${NC} GNU make not found\n"
    printf "${RED}error:${NC} PRIDE-PPPAR installation failed\n"; exit
fi

# Compilation & Installation
if [[ "${HOME}" == /root* ]]; then
    printf "${RED}error:${NC} unable to install PRIDE-PPPAR in /root \n"
    printf "${RED}error:${NC} PRIDE-PPPAR installation failed\n"; exit
fi

install_dir=${HOME}/.PRIDE_PPPAR_BIN
rm -rf "$install_dir"
cd src && make clean && make && make install \
    && cd .. \
    && mkdir -p $install_dir \
    && cp -f ./bin/* $install_dir \
    && chmod 755 ./scripts/*

# Check python
if which python3 > /dev/null 2>&1; then
    sed -i 's/\#\!\/usr\/bin\/env python$/\#\!\/usr\/bin\/env python3/g' ./scripts/*.py
elif which python > /dev/null 2>&1; then
    sed -i 's/\#\!\/usr\/bin\/env python3$/\#\!\/usr\/bin\/env python/g' ./scripts/*.py
else
    printf "${YELLOW}warning:${NC} no python or python3, please install one of them before using py scripts\n"
fi

cp -f ./scripts/*py ./scripts/*sh $install_dir
mv -f $install_dir/pdp3.sh $install_dir/pdp3
if [ $? -eq 0 ]; then
    if [ "$SYS" == "Darwin" ]; then
        grep "^export PATH=$install_dir:\$PATH" ${HOME}/.bash_profile > /dev/null 2>&1
        [ $? -ne 0 ] && echo "export PATH=$install_dir:\$PATH" >> ${HOME}/.bash_profile
    else
        grep "^export PATH=$install_dir:\$PATH" ${HOME}/.bashrc > /dev/null 2>&1
        [ $? -ne 0 ] && echo "export PATH=$install_dir:\$PATH" >> ${HOME}/.bashrc
    fi
fi

# Output
ls ${install_dir}/lsq > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo -e "\033[1;31m" && cat ./doc/logo && echo -e "$NC"
    printf "${BLUE}::${NC} PRIDE-PPPAR (v3.0) installation successfully completed!\n"
    printf "${BLUE}::${NC} executable binaries are copy to $install_dir\n"
    printf "${BLUE}::${NC} $install_dir added to PATH\n"
else
    printf "${RED}error:${NC} PRIDE-PPPAR installation failed!\n"
    exit 1
fi

chmod 755 $install_dir/*

work_dir=$(pwd)
cp ./table/config_template $install_dir
if [ "$SYS" == "Darwin" ]; then
    sed -i '' "/Table directory/s#Default#$work_dir/table/#" $install_dir/config_template
else
    sed -i    "/Table directory/s#Default#$work_dir/table/#" $install_dir/config_template
fi

chmod 444 $install_dir/config_template
chmod 755 ./example/*.sh

# Select download sites
printf "\n"
read -p $'Would you like to download GNSS products from IGN (\e[31mrecommended for users outside China\e[0m) [Y/N]: ' test
if [ ${#test} -ge 1 ] && ( [ ${test:0:1} == "y" ] || [ ${test:0:1} == "Y" ] ); then
    touch $install_dir/.ign_priority
fi

# Test examples
printf "\n"
read -p $'Run tests or not (\e[31mstrongly recommended for the first installation !!!\e[0m) [Y/N]: ' test
if [ ${#test} -ge 1 ] && ( [ ${test:0:1} == "y" ] || [ ${test:0:1} == "Y" ] ); then
    cd example
    if [ "$SYS" == "Darwin" ]; then
        source ~/.bash_profile
    else
        source ~/.bashrc
    fi
    bash test.sh
    cd ..
fi
