#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: GNSS PPP&PPP-AR data processing with PRIDE PPP-AR 2             ##
##                                                                           ##
##  AUTHOR : PRIDE LAB      pride@whu.edu.cn                                 ##
##                                                                           ##
##  VERSION: ver 2.1                                                         ##
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


######################################################################
##                        Message Colors                            ##
######################################################################
NC='\033[0m'
RED='\033[0;31m'
BLUE='\033[1;34m'
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

MSGERR="${RED}error:$NC"
MSGWAR="${YELLOW}warning:$NC"
MSGINF="${BLUE}::$NC"
MSGSTA="${BLUE}===>$NC"

USECACHE=YES    # YES/NO (uppercase!)
DEBUG=NO

######################################################################
##                     Funciton definations                         ##
######################################################################
main()
{
    CheckCmdArgs "$@" || exit 1
    CheckExecutables || exit 1

    local ctrl_file="$1"
    local ymd_start=(${2:0:4} ${2:4:2} ${2:6:2})
    local ymd_end=(${3:0:4} ${3:4:2} ${3:6:2})
    local site="${4:0:4}"
    local mode="${5:0:1}"
    local interval="$6"
    local AR="${7:0:1}"

    ctrl_file=$(readlink -f "$ctrl_file") # convert to absolute path

    # Output processing infomation
    echo -e "$MSGINF Processing date range: ${ymd_start[*]}  <==>  ${ymd_end[*]}"
    echo -e "$MSGINF Site name: ${site}"
    echo -e "$MSGINF Positioning mode: ${mode}"
    echo -e "$MSGINF Control file: ${ctrl_file}"
    echo -e "$MSGINF AR switch: ${AR}"

    # Processing day-by-day
    readonly local mjd_start=$(ymd2mjd ${ymd_start[*]})
    readonly local mjd_end=$(ymd2mjd ${ymd_end[*]})
    local work_dir=$(pwd) ydoy
    for mjd in $(seq $mjd_start $mjd_end)
    do
        cd "${work_dir}"
        ydoy=($(mjd2ydoy ${mjd}))
        mkdir -p ./${ydoy[0]}/${ydoy[1]}
        cd ./${ydoy[0]}/${ydoy[1]}
        if [ $? -eq 0 ]; then
            ProcessSingleDay $mjd "$ctrl_file" ${AR} ${site,,} ${mode^^} ${interval} || echo -e "$MSGWAR ${ydoy[*]} processing failed"
        else
            echo -e "$MSGERR no such directory: ${work_dir}/$year/$doy"
            echo -e "$MSGWAR skip processing: $year $doy"
        fi
    done
}

CheckCmdArgs() { # purpose: chech whether command line arguments are right
                 # usage  : CheckCmdArgs "$@"
    if [ $# -ne 7 ]; then
        PRIDE_PPPAR_Help
        return 1
    else
        local ctrl_file="$1"
        local ymd_start="$2"
        local ymd_end="$3"
        local site="$4"
        local mode="$5"
        local interval="$6"
        local AR="$7"
        if [ ! -e "$ctrl_file" ]; then
            echo -e "$MSGERR no control file: $ctrl_file"
            return 1
        elif [ ${#ymd_start} -ne 8 ]; then
            echo -e "$MSGERR bad start_date format: $ymd_start"
            echo -e "$MSGINF start_date format: YYYYMMDD"
            return 1
        elif [ ${#ymd_end} -ne 8 ]; then
            echo -e "$MSGERR bad end_date format: $ymd_start"
            echo -e "$MSGINF end_date format: YYYYMMDD"
            return 1
        elif [ ${#site} -ne 4 ]; then
            echo -e "$MSGERR bad site name format: $site"
            echo -e "$MSGINF site name format: NNNN(nnnn)"
            return 1
        elif [ ${mode,,} != s -a ${mode,,} != k -a ${mode,,} != f ]; then
            echo -e "$MSGERR unknown positioning mode: $mode"
            echo -e "$MSGINF positioning mode: S(s)/K(k)/F(f)"
            return 1
        elif [ $(echo "${interval} < 0.02"|bc) -eq 1 -o  $(echo "${interval} > 30"|bc) -eq 1 ]; then
            echo -e "$MSGERR too small/big interval: $interval"
            echo -e "$MSGINF note: 0.02<=interval<=30"
            return 1
        elif [ ${AR,,} != y -a ${AR,,} != n ]; then
            echo -e "$MSGERR unknown AR option: $AR"
            echo -e "$MSGINF AR option: Y(y)/N(n)"
            return 1
        fi
    fi
}

CheckExecutables() { # purpose: check whether all needed executables are callable
                     # usage  : CheckExecutables
    echo -e "$MSGSTA CheckExecutables..."
    if [ which lsq > /dev/null 2>&1 ]; then
        echo -e "$MSGERR lsq not found" && return 1
    fi
    if [ which arsig > /dev/null 2>&1 ]; then
        echo -e "$MSGERR arsig not found" && return 1
    fi
    if [ which tedit > /dev/null 2>&1 ]; then
        echo -e "$MSGERR tedit not found" && return 1
    fi
    if [ which redig > /dev/null 2>&1 ]; then
        echo -e "$MSGERR redig not found" && return 1
    fi
    if [ which get_ctrl > /dev/null 2>&1 ]; then
        echo -e "$MSGERR get_ctrl not found" && return 1
    fi
    if [ which spp > /dev/null 2>&1 ]; then
        echo -e "$MSGERR spp not found" && return 1
    fi
    if [ awk --help > /dev/null 2>&1 ]; then
        echo -e "$MSGERR awk not found" && return 1
    fi
    if [ which wget > /dev/null 2>&1 ]; then
        echo -e "$MSGERR wget not found" && return 1
    fi
    if [ which readlink > /dev/null 2>&1 ]; then
        echo -e "$MSGERR readlink not found" && return 1
    fi
    echo -e "$MSGSTA CheckExecutables done"
}

PRIDE_PPPAR_Help() { # purpose: print usage for PRIDE_PPPAR
                     # usage  : PRIDE_PPPAR_Help
    echo " -----------------------------------------------------------------------"
    echo "  Purpose  :    GNSS PPP&PPP-AR data processing with PRIDE PPP-AR 2"
    echo "  Usage    :    pdp3 ctrl_file start_date end_date site mode(S/K/F) interval AR(Y/N)"
    echo "                   -- ctrl_file  : configuration file of PRIDE PPP-AR 2"
    echo "                   -- start_date : start date for processing, format: YYYYMMDD"
    echo "                   -- end_date   : end date for procesing, format: YYYYMMDD"
    echo "                   -- site       : site name, format: NNNN"
    echo "                   -- mode(S/K/F): positioning mode"
    echo "                   -- interval   : 0.02<=interval<=30"
    echo "                   -- AR(Y/N)    : Switch of Ambiguity Resolution"
    echo "  Example  :    pdp3 config 20200101 20200101 abmf S 30 Y"
    echo "  Copyright:    GNSS Research Center, Wuhan University, 2021"
    echo " -----------------------------------------------------------------------"
}

ProcessSingleDay() { # purpose: process data of single day
                     # usage  : ProcessSingleDay mjd ctrl_file AR(y/n) site mode(s/k/f)
    local mjd=$1
    local ctrl_file="$2"
    local AR=$3
    local site=$4
    local mode=$5
    local interval=$6

    local ydoy=($(mjd2ydoy $mjd))
    local year=${ydoy[0]}
    local doy=${ydoy[1]}
    local ymd=($(ydoy2ymd ${ydoy[*]}))
    local mon=${ymd[1]}
    local day=${ymd[2]}
    #unset ydoy ymd

    echo -e "$MSGSTA ProcessSingleDay $year $doy..."

    # Clean the directory
    rm -f sit.xyz igserp config\.*
    local tp types=(rck ztd htg amb res stt con neq att fcb orb sck)
    for tp in ${types[*]}
    do
        rm -f ${tp}_${year}${doy}
    done
    types=(log pos kin)
    for tp in ${types[*]}
    do
        rm -f ${tp}_${year}${doy}_${site}
    done

    # Create config for current day
    config=$(basename `mktemp -u`)
    config=${config/tmp/config}    # configuration file
    cp -f "$ctrl_file" $config || return 1

    sed -i -e "/^Session time/s/-YYYY-/$year/g" \
           -e "/^Session time/s/-MM-/$mon/g" \
           -e "/^Session time/s/-DD-/$day/g" "$config"
    sed -i "/Rinex directory/s/-YEAR-/$year/g; s/-DOY-/$doy/g" "$config"

    # Prepare tables
    local table_dir=$(get_ctrl "$config" "Table directory")
    local product_dir=$(get_ctrl "$config" "Product directory")
    CopyTables "$table_dir" $mjd || return 1
    # leap.sec
    local leapsec_ftp=0 leapsec_exi=0
    if [ -f leap.sec ]; then
        sed -n '1 p' leap.sec | grep "*" > /dev/null 2>&1
        leapsec_ftp=$?
    else
        leapsec_exi=$?
    fi
    if [ "${leapsec_ftp}" != 0 -o "${leapsec_exi}" != 0 ]; then
        local leapsec="leap.sec"
        local leapsec_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/table/${leapsec}"
        rm -f $leapsec
        WgetDownload "$leapsec_url"
        if [ $? -ne 0 ]; then
            cp -f ${table_dir}/${leapsec} ./
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR ProcessSingleDay: no leap.sec"
                return 1
            fi
        fi
    fi
    
    # Download rinexnav file
    local rinex_dir=$(get_ctrl "$config" "Rinex directory")
    local rinexnav="BRDC00IGS_R_${year}${doy}0000_01D_MN.rnx"
    local urlnav="ftp://igs.gnsswhu.cn/pub/gps/data/daily/${year}/$doy/${year:2:2}p/${rinexnav}.gz"
    if [ ! -f "${rinex_dir}/brdm${doy}0.${year:2:2}p" ]; then
        echo -e "$MSGSTA Downloading ${rinexnav}..."
        WgetDownload "$urlnav"
        if [ $? -eq 0 ]; then
            gunzip -f ${rinexnav}.gz && mv ${rinexnav} ${rinex_dir}/brdm${doy}0.${year:2:2}p
            [ -f ${rinex_dir}/brdm${doy}0.${year:2:2}p ] && echo -e "$MSGSTA Downloading ${rinexnav} done" || return 1
        else
            rinexnav="BRDC00IGN_R_${year}${doy}0000_01D_MN.rnx"
            urlnav="ftp://igs.ign.fr/pub/igs/data/${year}/$doy/${rinexnav}.gz"
            WgetDownload "$urlnav"
            if [ $? -eq 0 ]; then
                gunzip -f ${rinexnav}.gz && mv ${rinexnav} ${rinex_dir}/brdm${doy}0.${year:2:2}p
                [ -f ${rinex_dir}/brdm${doy}0.${year:2:2}p ] && echo -e "$MSGSTA Downloading ${rinexnav} done" || return 1
            else
                echo -e "$MSGERR download rinexnav failed: ${rinex_dir}/${rinexnav}"
                return 1
            fi
        fi
    fi

    # RINEX-OBS check
    rinexnav="${rinex_dir}/brdm${doy}0.${year:2:2}p"
    local rinexobs="${rinex_dir}/${site}${doy}0.${year:2:2}o"
    if [ ! -f "$rinexobs" ]; then
        echo -e "$MSGWAR ${rinexobs} doesn't exist"
        rinexobs="${rinex_dir}/${site^^}00[A-Z][A-Z][A-Z]_R_${year}${doy}0000_01D_30S_MO.rnx"
        local n_rinexobs=$(ls ${rinexobs} 2> /dev/null | wc -l)
        if [ "$n_rinexobs" -ne 1 ]; then
            echo -e "$MSGWAR ${rinex_dir}/${site^^}00XXX_R_${year}${doy}0000_01D_30S_MO.rnx doesn't exist" && return 1
        fi
    fi

    # Compute a priori positions
    echo -e "$MSGSTA Prepare initial position ${site}..."
    ComputeInitialPos "$rinexobs" "$rinexnav" "$interval" "$mode"
    awk -v sit=$site '{if(NR==1){printf(" %s%16.4f%16.4f%16.4f\n",sit,$3,$4,$5)}}' tmp_ComputeInitialPos > sit.xyz
    #local obs_session=($(awk '{if(NR==2){printf("%04d %02d %02d %02d %02d %02d %7d\n",$3,$4,$5,$6,$7,$8,$9)}}' tmp_ComputeInitialPos))
    local obs_session=($(awk '{if(NR==2){print $3,$4,$5,$6,$7,$8,$9}}' tmp_ComputeInitialPos))
    rm -f tmp_ComputeInitialPos
    if [ ${mode} == "F" ]; then
        local initial_pos=($(snxsit.sh $site ${ymd[*]}))
        if [ ${#initial_pos[@]} -ne 6 ]; then
            echo -e "$MSGERR ProcessSingleSite: no position or sigma: $site"
            return 1
        fi
        printf " %s%16.4f%16.4f%16.4f%10.6f%10.6f%10.6f\n" ${site} ${initial_pos[*]} > sit.xyz
    fi
    echo -e "$MSGSTA Prepare initial position ${site} done"

    # Fill in session times
    if [ ${#obs_session[@]} -ne 7 ]; then
        echo -e "$MSGERR ProcessSingleDay: no session time"
        return 1
    else
        sed -i -e "/^Session time/s/-HH-/${obs_session[3]}/g" \
               -e "/^Session time/s/-MI-/${obs_session[4]}/g" \
               -e "/^Session time/s/-SS-/${obs_session[5]}/g" \
               -e "/^Session time/s/-SE-/${obs_session[6]}/g" "$config"
    fi

    # Prepare products
    PrepareProducts $mjd "$product_dir" ${config} ${AR}
    if [ $? -ne 0 ]; then
        echo -e "$MSGERR PrepareProducts failed"
        return 1
    fi

    # process single day
    sites=($(awk '/^ [0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z] [KSF]/ {print $1}' "$config"))
    [ ${#sites[@]} -eq 0 ] && echo -e "$MSGWAR: no option line to be processed" && return 1
    [ ${#sites[@]} -gt 1 ] && echo -e "$MSGWAR: more than one option line to be processed" && return 1
    sed -i "s/^ \w\w\w\w [KSF]/ ${site} ${mode}/" ${config}
    sed -i "/Interval/s/[0-9]\.*[0-9]*/$interval/" ${config}
    ProcessSingleSite "$site" $year $doy "${config}" $AR
    if [ $? -ne 0 ]; then
        echo -e "$MSGWAR ProcessSingleDay: skip processing $year $doy $site"
        if [ ${DEBUG} == "NO" ]; then
            types=(rck ztd htg amb res stt con neq)
            for tp in ${types[*]}
            do
                rm -f ${tp}_${year}${doy}
            done
            types=(log pos kin)
            for tp in ${types[*]}
            do
                rm -f ${tp}_${year}${doy}_${site}
            done
        fi
    else
        echo -e "$MSGSTA ProcessSingleDay $year $doy done"
    fi
    #rm -f ${config}
}

ProcessSingleSite() { # purpose: process data of single site
                      # usage  : ProcessSingleSite site year doy config AR(Y/N)
    local site=$1
    local year=$2
    local doy=$3
    local config=$4
    local AR=${5:0:1}
    local position_mode=$(grep "^ $site [KSF]" "$config" | awk '{print $2}') # Static/Kinematic
    local cutoff_elev=$( grep "^ $site [KSF]" "$config" | awk '{print $5}') # int, degree

    echo -e "$MSGSTA ProcessSingleSite $year $doy $site..."

    # Create kin file for K mode for spp
    local editing=$(get_ctrl "$config" "Strict editing")
    if [ "$editing" == "YES" ]; then
        local editing_mode="YES"
    elif [ "$editing" == "NO" ]; then
        local editing_mode="NO"
    else
        echo -e "$MSGERR ProcessSingleSite: unknown editing mode: $editing"
        return 1
    fi

    # Data preprocess
    echo -e "$MSGSTA Data pre-processing..."
    local interval=$(get_ctrl "$config" "Interval")
    local session=$(grep "^Session time" "${config}" | awk '{print $10}')
    local hms=($(grep "^Session time" "${config}" | awk '{print $7,$8,$9}'))
    local rhd_file="log_${year}${doy}_${site}"
    local ymd=($(ydoy2ymd $year $doy))
    xyz=($(awk -v sit=$site '{if($1==sit){print $2,$3,$4}}' sit.xyz))
    local cmd=""
    if [ "$position_mode" == S -o "$position_mode" == F ]; then
        cmd="tedit ${rinexobs} -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
            -xyz ${xyz[*]} -short 1200 -lc_check only -rhd ${rhd_file} -pc_check 300 \
            -elev ${cutoff_elev} -rnxn ${rinexnav}"
        local mjd=$(ymd2mjd ${ymd[*]})
        if [ $mjd -le 51666 ]; then
            cmd="tedit ${rinexobs} -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
                -xyz ${xyz[*]} -short 1200 -lc_check no -rhd ${rhd_file} -pc_check 0 \
                -elev ${cutoff_elev} -rnxn ${rinexnav}"
        fi
    elif [ "$position_mode" == K ]; then
        cmd="tedit ${rinexobs} -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
              -xyz kin_${year}${doy}_${site} -short 120 -lc_check no \
             -elev ${cutoff_elev} -rhd ${rhd_file} -rnxn ${rinexnav}"
        local mjd=$(ymd2mjd ${ymd[*]})
        if [ $mjd -le 51666 ]; then
            cmd="tedit ${rinexobs} -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
                 -xyz kin_${year}${doy}${site} -short 120 -lc_check no \
                 -pc_check 0 -elev ${cutoff_elev} -rhd ${rhd_file} -rnxn ${rinexnav}"
        fi
    else
        echo -e "$MSGERR ProcessSingleSite: unknown position mode: $site $position_mode"
        return 1
    fi
    cmd=$(tr -s " " <<< "$cmd")
    ExecuteWithoutOutput "$cmd" || return 1
    echo -e "$MSGSTA Data pre-processing done"

    # Data clean (iteration)
    echo -e "$MSGSTA Data cleaning..."
    if [ "$editing_mode" == YES ]; then
        local short=$(echo $interval | awk '{printf("%.0f\n", 600/$1)}')
        local jumps=(400 200 100 50)
        local jump_end=50
    else
        local short=$(echo $interval | awk '{printf("%.0f\n", 120/$1)}')
        local jumps=(400 200 100)
        local jump_end=100
    fi
    for jump in ${jumps[*]}
    do
        cmd="lsq ${config}"
        ExecuteWithoutOutput "$cmd" || return 1
        cmd="redig res_${year}${doy} -jmp $jump -sht $short"
        local time=`date +'%Y-%m-%d %H:%M:%S'`
        $cmd > tempout 2>&1
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} executed ok"
        else
            echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} executed failed"
            return 1
        fi
        awk '/%%%\+RMS OF RESIDUALS---PHASE\(MM\)/,/%%%\-RMS OF RESIDUALS---PHASE\(MM\)/{print}' tempout
        new_rem=`awk '/NEWLY REMOVED:/{print $3}' tempout`
        awk '/NEWLY REMOVED:/{printf "\033[1;34mNewly removed observations\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
        new_amb=`awk '/NEWLY AMBIGUT:/{print $3}' tempout`
        awk '/NEWLY AMBIGUT:/{printf "\033[1;34mNewly inserted ambiguities\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
    done
    njump=0
    while [ $new_rem != 0 -o $new_amb != 0 ]
    do
        [ $njump -gt 100 ] && break
        ((njump=njump+1))
        cmd="lsq ${config}"
        ExecuteWithoutOutput "$cmd" || return 1
        cmd="redig res_${year}${doy} -jmp $jump_end -sht $short"
        local time=`date +'%Y-%m-%d %H:%M:%S'`
        $cmd > tempout 2>&1
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} executed ok"
        else
            echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} executed failed"
            return 1
        fi
        new_rem=`awk '/NEWLY REMOVED:/{print $3}' tempout`
        new_amb=`awk '/NEWLY AMBIGUT:/{print $3}' tempout`
        if [ $new_rem == '' -o $new_amb == '' ]; then
            echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} executed failed"
            return 1
        fi
        awk '/%%%\+RMS OF RESIDUALS---PHASE\(MM\)/,/%%%\-RMS OF RESIDUALS---PHASE\(MM\)/{print}' tempout
        awk '/NEWLY REMOVED:/{printf "\033[1;34mNewly removed observations\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
        awk '/NEWLY AMBIGUT:/{printf "\033[1;34mNewly inserted ambiguities\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
    done
    rm -f tempout
    echo -e "$MSGSTA Data cleaning done"

    # Ambiguity fixing
    if [ $AR == Y -o $AR == y ]; then
        cmd="arsig ${config}"
        Execute "$cmd" || return 1
        cmd="lsq ${config}"
        Execute "$cmd" || return 1
    fi
    echo -e "$MSGSTA Final processing done"

    # Rename result files
    local tp types=(rck ztd htg amb res stt con) fn
    for tp in ${types[*]}
    do
        fn=${tp}_${year}${doy}
        [ -f ${fn} ] && mv -f ${fn} ${fn}_${site}
    done

    echo -e "$MSGSTA ProcessSingleSite $year $doy $site done"
}

ComputeInitialPos() { # purpose: compute intial postion with spp
                      # usage  : ComputeInitialPos rinexobs rinexnav
    local rinexobs="$1"
    local rinexnav="$2"
    local interval="$3"
    local mode="$4"
    #local site=$(basename "$rinexobs"); site=${site:0:4}

    local cmd=""
    if [ "$mode" == "K" ]; then
        cmd="spp -trop "saas" -ti ${interval} -o kin_${year}${doy}_${site} ${rinexobs} ${rinexnav}"
    else
        cmd="spp -trop "saas" -ti ${interval} ${rinexobs} ${rinexnav}"
    fi
    Execute "$cmd" tmp_ComputeInitialPos || return 1
}

CopyTables() { # purpose: copy PRIDE-PPPAR needed tables to working directory
               # usage  : CopyTables table_dir mjd
    echo -e "$MSGSTA CopyTables..."
    local table_dir="$1"
    local mjd=$2
    local tables=(file_name oceanload orography_ell orography_ell_1x1 gpt3_1.grd)
    for table in ${tables[*]}
    do
        if [ ! -f "$table_dir/$table" ]; then
             echo -e "$MSGERR CopyTables: no such file: $table_dir/$table"
             return 1
        fi
        ln -sf "$table_dir/$table" .
    done

    echo -e "$MSGSTA CopyTables done"
}

PrepareProducts() { # purpose: prepare PRIDE-PPPAR needed products in working directory
                    # usage  : PrepareProducts mjd products_dir AR(y/n)
    local mjd_mid=$1
    local products_dir="$2"
    local config="$3"
    local AR=${4:0:1}
    [ -d $products_dir ] || mkdir -p "$products_dir"

    echo -e "$MSGSTA PrepareProducts..."

    local ydoy=($(mjd2ydoy $mjd_mid))
    local ymd=($(ydoy2ymd ${ydoy[*]}))
    local wkdow=($(mjd2wkdow $mjd_mid))
    local year=${ydoy[0]}
    local doy=${ydoy[1]}

    local custom_pro_clk=$(get_ctrl "$config" "Satellite clock")
    if [ "$custom_pro_clk" != Default ]; then        
        local clk="${custom_pro_clk}"
        local clk_url="${clk}"
        local argnum=`get_ctrl "$config" "Satellite clock" | wc -w`
        if [ $argnum -eq 1 ]; then
            CopyOrDownloadProduct "$products_dir/$clk" "$clk_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: no such file: $clk"
                return 1
            fi          
        elif [ $argnum -gt 1 ]; then    
            MergeProducts "$products_dir" "$clk" "mersck_${year}${ydoy[1]}"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: $clk merge failed"
                return 1
            fi
            clk="mersck_${year}${ydoy[1]}"
        fi        
    else
        local clk="WUM0MGXRAP_${year}${ydoy[1]}0000_01D_30S_CLK.CLK.gz"
        local clk_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/clock/${clk}"
        CopyOrDownloadProduct "$products_dir/$clk" "$clk_url"
        if [ $? -ne 0 ]; then
            echo -e "$MSGERR PrepareProducts: $clk download failed"
            return 1
        else
            [ -f ${clk} ] && gunzip -f ${clk}
        fi
    fi
    [ "$custom_pro_clk" == Default ] && sed -i "/Satellite clock/s/Default/${clk}/" $config

    local custom_pro_fcb=$(get_ctrl "$config" "Code/phase bias")
    if [ "$custom_pro_fcb" != Default ]; then
        local fcb="${custom_pro_fcb}"
        local fcb_url="${fcb}"
        local argnum=`get_ctrl "$config" "Code/phase bias" | wc -w`
        if [ $argnum -eq 1 ]; then            
            CopyOrDownloadProduct "$products_dir/$fcb" "$fcb_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGWAR PrepareProducts: no such file: $fcb"
                [ $AR == y -o $AR == Y ] && echo -e "$MSGERR no phase bias product: $fcb" && return 1
            fi
        elif [ $argnum -gt 1 ]; then           
            MergeProducts "$products_dir" "$fcb" "merfcb_${year}${ydoy[1]}"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: $fcb merge failed"
                return 1
            fi
            fcb="merfcb_${year}${ydoy[1]}"
        fi
    else
        local fcb="WUM0MGXRAP_${year}${ydoy[1]}0000_01D_01D_ABS.BIA.gz"
        local fcb_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/bias/${fcb}"
        CopyOrDownloadProduct "$products_dir/$fcb" "$fcb_url"
        if [ $? -ne 0 ]; then
            echo -e "$MSGWAR PrepareProducts: $fcb download failed"
            [ $AR == y -o $AR == Y ] && echo -e "$MSGERR no phase bias product: $fcb" && return 1
        else
            [ -f ${fcb} ] && gunzip -f ${fcb}
        fi
    fi
    [ "$custom_pro_fcb" == Default ] && sed -i "/Code\/phase bias/s/Default/${fcb}/" $config
    
    # IGS ATX
    echo -e "$MSGINF Prepare IGS ATX product: $abs_atx ..."
    local abs_atx=null
    if [ "$custom_pro_clk" != Default ]; then
        abs_atx="igs14.atx"
    else
        abs_atx="$(grep "SYS / PCVS APPLIED$" ${clk%.gz}|head -1|cut -c21-34|tr 'A-Z' 'a-z')" 
    fi       
    local atx_url="https://files.igs.org/pub/station/general/pcv_archive/${abs_atx}"
    CopyOrDownloadProduct "$table_dir/$abs_atx" "$atx_url"
    if [ $? -ne 0 ]; then
        atx_url="https://files.igs.org/pub/station/general/${abs_atx}"
        CopyOrDownloadProduct "$table_dir/$abs_atx" "$atx_url"
        if [ $? -ne 0 ]; then
            echo -e "$MSGERR PrepareProducts: no such file: $table_dir/$abs_atx"
            return 1
        fi
    fi
    ln -sf $table_dir/$abs_atx . && mv $abs_atx abs_igs.atx
    echo -e "$MSGINF Prepare IGS ATX product: $abs_atx done"

    local custom_pro_sp3=$(get_ctrl "$config" "Satellite orbit")
    if [ "$custom_pro_sp3" != Default ]; then
        local sp3="${custom_pro_sp3}"
        local sp3_url="${sp3}"
        local argnum=`get_ctrl "$config" "Satellite orbit" | wc -w`
        if [ $argnum -eq 1 ]; then           
            CopyOrDownloadProduct "$products_dir/$sp3" "$sp3_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: no such file: $sp3"
                return 1
            fi
        elif [ $argnum -gt 1 ]; then
            MergeProducts "$products_dir" "$sp3" "mersp3_${year}${ydoy[1]}"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: $sp3 merge failed"
                return 1
            fi
            sp3="mersp3_${year}${ydoy[1]}"
        fi
    else
        local sp3="WUM0MGXRAP_${year}${ydoy[1]}0000_01D_01M_ORB.SP3.gz"
        local sp3_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/orbit/${sp3}"
        CopyOrDownloadProduct "$products_dir/$sp3" "$sp3_url"
        if [ $? -ne 0 ]; then
            echo -e "$MSGERR PrepareProducts: $sp3 download failed"
            return 1
        else
            [ -f ${sp3} ] && gunzip -f ${sp3}
        fi
    fi
    [ "$custom_pro_sp3" == Default ] && sed -i "/Satellite orbit/s/Default/${sp3}/" $config

    local custom_pro_erp=$(get_ctrl "$config" "ERP")
    if [ "$custom_pro_erp" != Default ]; then
        local erp="${custom_pro_erp}"
        local erp_url="${erp}"
        local argnum=`get_ctrl "$config" "ERP" | wc -w`
        if [ $argnum -eq 1 ]; then
            CopyOrDownloadProduct "$products_dir/$erp" "$erp_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: no such file: $erp"
                return 1
            fi
        elif [ $argnum -gt 1 ]; then            
            MergeProducts "$products_dir" "$erp" "mererp_${year}${ydoy[1]}"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: $erp merge failed"
                return 1
            fi
            erp="mererp_${year}${ydoy[1]}"
        fi
    else
        local erp="WUM0MGXRAP_${year}${ydoy[1]}0000_01D_01D_ERP.ERP.gz"
        local erp_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/orbit/${erp}"
        CopyOrDownloadProduct "$products_dir/$erp" "$erp_url"
        if [ $? -ne 0 ]; then
            echo -e "$MSGERR PrepareProducts: $erp download failed"
            return 1
        else
            [ -f ${erp} ] && gunzip -f ${erp}
        fi
    fi
    [ "$custom_pro_erp" == Default ] && sed -i "/ERP/s/Default/${erp}/" $config

    local custom_pro_att=$(get_ctrl "$config" "Quaternions")
    if [ "$custom_pro_att" != Default ]; then
        local att="${custom_pro_att}"
        local att_url="${att}"
        local argnum=`get_ctrl "$config" "Quaternions" | wc -w`
        if [ $argnum -eq 1 ]; then
            CopyOrDownloadProduct "$products_dir/$att" "$att_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGWAR PrepareProducts: no such file: $att"
            fi
        elif [ $argnum -gt 1 ]; then           
            MergeProducts "$products_dir" "$att" "meratt_${year}${ydoy[1]}"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: $att merge failed"
                return 1
            fi
            att="meratt_${year}${ydoy[1]}"
        fi
    else
        local att="WUM0MGXRAP_${year}${ydoy[1]}0000_01D_30S_ATT.OBX.gz"
        local att_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/orbit/${att}"
        CopyOrDownloadProduct "$products_dir/$att" "$att_url"
        if [ $? -ne 0 ]; then
            echo -e "$MSGWAR PrepareProducts: $att download failed"
        else
            [ -f ${att} ] && gunzip -f ${att}
        fi
    fi
    [ "$custom_pro_att" == Default ] && sed -i "/Quaternions/s/Default/${att}/" $config

    local ion_file=$(get_ctrl "$config" "Iono 2nd")
    if [ ${ion_file} == "YES" ]; then
        echo -e "$MSGSTA Downloading High-order Ion Grid..."
        local ion ion_url
        ion="CODG${doy}0.${year:2:2}I.Z"
        ion_url="ftp://ftp.aiub.unibe.ch/CODE/${year}/${ion}"
        CopyOrDownloadProduct "$products_dir/$ion" "$ion_url" || return 1
        [ -f ${ion} ] && gunzip -f ${ion} && mv ${ion%.Z} tec_${ydoy[0]}${ydoy[1]}
        echo -e "$MSGSTA Downloading High-order Ion Grid done"
    fi
    
    grep '^ [0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z] .*VM1' ${config} 2>&1 > /dev/null
    if [ $? -eq 0 ]; then
        echo -e "$MSGSTA Downloading VMF1 GRID..."
        local vmf vmf_url hour
        # Current Day (for interpolation)
        for hour in `seq 0 6 18 | awk '{printf("%02d\n",$1)}'`
        do
            vmf="VMFG_${ymd[0]}${ymd[1]}${ymd[2]}.H${hour}"
            vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${ydoy[0]}/${vmf}"
            CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1
        done

        # Previous Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_mid-1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMFG_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H18"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1

        # Next Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_mid+1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMFG_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H00"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1

        cat VMFG_* > vmf_${ydoy[0]}${ydoy[1]} || return 1
        echo -e "$MSGSTA Downloading VMF1 GRID done"
    fi

    grep '^ [0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z] .*VM3' ${config} 2>&1 > /dev/null
    if [ $? -eq 0 ]; then
        echo -e "$MSGSTA Downloading VMF3 GRID..."
        local vmf vmf_url hour
        # Current Day (for interpolation)
        for hour in `seq 0 6 18 | awk '{printf("%02d\n",$1)}'`
        do
            vmf="VMF3_${ymd[0]}${ymd[1]}${ymd[2]}.H${hour}"
            vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/${ydoy[0]}/${vmf}"
            CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1
        done

        # Previous Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_mid-1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMF3_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H18"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1

        # Next Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_mid+1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMF3_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H00"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1

        cat VMF3_* > vmf_${ydoy[0]}${ydoy[1]} || return 1
        echo -e "$MSGSTA Downloading VMF3 GRID done"
    fi

    # rename products
    if [ "$custom_pro_clk" != Default ]; then
        mv ${clk} sck_${ydoy[0]}${ydoy[1]} || return 1
    else
        mv ${clk%.gz} sck_${ydoy[0]}${ydoy[1]} || return 1
    fi
    if [ "$custom_pro_fcb" != Default ]; then
        [ -e ${fcb} ] && mv ${fcb} fcb_${ydoy[0]}${ydoy[1]}
    else
        [ -e ${fcb%.gz} ] && mv ${fcb%.gz} fcb_${ydoy[0]}${ydoy[1]}
    fi
    if [ "$custom_pro_att" != Default ]; then
        [ -e ${att} ] && mv ${att} att_${ydoy[0]}${ydoy[1]}
    else
        [ -e ${att%.gz} ] && mv ${att%.gz} att_${ydoy[0]}${ydoy[1]}
    fi

    # Generate igserp
    if [ "$custom_pro_erp" != Default ]; then
        mv ${erp} igserp || return 1
    else
        mv ${erp%.gz} igserp || return 1
    fi

    # Generate binary sp3
    if [ "$custom_pro_sp3" != Default ]; then
        cmd="sp3orb ${sp3} -cfg ${config}"
        ExecuteWithoutOutput "${cmd}" || return 1
    else
        cmd="sp3orb ${sp3%.gz} -cfg ${config}"
        ExecuteWithoutOutput "${cmd}" || return 1
    fi  

    echo -e "$MSGSTA PrepareProducts done"
}

MergeProducts() {
    local dir="$1"
    local infile="$2"
    local outfile="$3"
    rm -f $outfile
    for f in $infile
    do
        cat "$dir/$f" >> $outfile
    done
    [ -f $outfile ] && return 0 || return 1
}

CopyOrDownloadProduct() { # purpose: copy or download a product
                          # usage  : CopyOrDownloadProduct file url
    local file="$1"
    local url="$2"
    if [ ${USECACHE} = "YES" -a -f "$file" ]; then
        cp -f "$file" .
    else
        WgetDownload "$url" || return 1
        cp -f $(basename "$url") "$file"
        return 0
    fi
}

WgetDownload() { # purpose: download a file with wget
                 # usage  : WgetDownload url
    local url="$1"
    local args="-q -nv -nc -c -t 3 --connect-timeout=10 --read-timeout=60 --show-progress"
    cmd="wget ${args} ${url}"
    $cmd
    [ -e $(basename "${url}") ] && return 0 || return 1
}

LastYearMonth() { # purpose: get last year-month
                  # usage  : LastYearMonth year month
    local year=$1
    local mon=$((10#$2))
    [ $((mon-1)) -lt 1  ] && mon=12 && year=$((year-1)) || mon=$((mon-1))
    printf "%4d %02d\n" $year $mon
}

UncompressFile() { # purpose: uncompress a file automatically
                   # usage  : UncompressFile file del(Y/N)
    local file="$1"
    local del=$2
    # file $file
}

Execute() {
    local cmd="$1"
    if [ $# -gt 1 ]; then
        local outp="$2"
    fi
    # echo $cmd
    time=`date +'%Y-%m-%d %H:%M:%S'`
    if [ $# -gt 1 ]; then
        $cmd > $outp
    else
        $cmd
    fi
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} executed ok" # | tee -a $log
        return 0
    else
        echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} executed failed"  #| tee -a $log
        return 1
    fi
}

ExecuteWithoutOutput() {
    local cmd="$1"
    # echo $cmd
    time=`date +'%Y-%m-%d %H:%M:%S'`
    $cmd > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} executed ok" # | tee -a $log
        return 0
    else
        echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} executed failed"  #| tee -a $log
        echo -e "$MSGINF Here is the output:\n"
        $cmd
        return 1
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

mjd2ydoy()
{
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

mjd2wkdow()
{
    local mjd=$1
    local mjd0=44243
    local difmjd=$(($mjd-$mjd0-1))
    local week=$(($difmjd/7))
    local dow=$(($difmjd%7))
    printf "%04d %d\n" $week $dow
}

ydoy2ymd()
{
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


######################################################################
##                               Entry                              ##
######################################################################
main "$@"

