#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: GNSS PPP&PPP-AR data processing with PRIDE PPP-AR 3             ##
##                                                                           ##
##  AUTHOR : PRIDE LAB      pride@whu.edu.cn                                 ##
##                                                                           ##
##  VERSION: ver 3.0                                                         ##
##                                                                           ##
##  DATE   : Mar-08, 2024                                                    ##
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

######################################################################
##                          Message Color                           ##
######################################################################

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

######################################################################
##                          Basic Setting                           ##
######################################################################

shopt -s extglob                        # Enable Extendded Globbing
readonly LC_NUMERIC="en_US.UTF-8"       # Specify period decimal point
readonly OS="$(uname)"                  # Operation System

readonly DEBUG=YES                      # YES/NO (uppercase!)
readonly OFFLINE=NO                     # OFFLINE=YES will overwrite USECACHE=NO
readonly USECACHE=YES
readonly USERTS=YES

readonly SCRIPT_NAME="pdp3"
readonly VERSION_NUM="3.0"

######################################################################
##                     System-specific Command                      ##
######################################################################

sedi() {
    local cmd="$1"
    local fil="$2"
    if [ "$OS" == "Darwin" ]; then
        sed -i '' "$cmd" "$fil"
    else
        sed -i "$cmd" "$fil"
    fi
}

rdlk() {
    local arg="$1"
    if [ "$OS" == "Darwin" ]; then
        [[ $1 = /* ]] && echo "$arg" || echo "$(pwd)/${arg#./}"
    else
        readlink -f "$arg"
    fi
}

######################################################################
##                       Funciton Defination                        ##
######################################################################

main() {
    args=$(ParseCmdArgs "$@") || exit 1
    CheckExecutables          || exit 1

    local rnxo_path=$(echo "$args" | sed -n 1p)       # absolute path
    local ctrl_path=$(echo "$args" | sed -n 2p)       # absolute path
    local ctrl_file=$(echo "$args" | sed -n 3p)       # temporary config
    local    date_s=$(echo "$args" | sed -n 4p)       # yyyy-mm-dd
    local    hour_s=$(echo "$args" | sed -n 5p)       # hh:mi:ss
    local    date_e=$(echo "$args" | sed -n 6p)       # yyyy-mm-dd
    local    hour_e=$(echo "$args" | sed -n 7p)       # hh:mi:ss
    local      mode=$(echo "$args" | sed -n 8p)       # S/P/K/F/L, upper case
    local        AR=$(echo "$args" | sed -n 9p)       # A/Y/N, upper case

    local interval=$(get_ctrl "$ctrl_file" "Interval")
    local freq_cmb=$(get_ctrl "$ctrl_file" "Frequency combination")
    local site=$(grep "^ .... [A-Z]" "$ctrl_file" | cut -c 2-5)

    local rinex_dir=$(dirname  "$rnxo_path")
    local rnxo_name=$(basename "$rnxo_path")

    # Output processing infomation
    echo -e "$MSGINF Processing time range: $date_s $hour_s <==> $date_e $hour_e"
    echo -e "$MSGINF Processing interval: $interval"
    echo -e "$MSGINF Site name: $site"
    echo -e "$MSGINF Positioning mode: $mode"
    echo -e "$MSGINF AR switch: $AR"
    echo -e "$MSGINF Frequency combination: $freq_cmb"
    echo -e "$MSGINF Configuration file: $ctrl_path"
    echo -e "$MSGINF RINEX observation file: $rnxo_path"

    
    if [ "$OS" == "Darwin" ]; then
        local doy_s=$(date -j -f "%Y-%m-%d" "$date_s" +"%j")
        local doy_e=$(date -j -f "%Y-%m-%d" "$date_e" +"%j")
    else
        local doy_s=$(date -d "$date_s" +"%j")
        local doy_e=$(date -d "$date_e" +"%j")
    fi
    local ymd_s=($(echo "$date_s" | tr '-' ' '))
    local ymd_e=($(echo "$date_e" | tr '-' ' '))
    local mjd_s=$(ymd2mjd ${ymd_s[*]})
    local mjd_e=$(ymd2mjd ${ymd_e[*]})

    local proj_dir=$(pwd)
    local mjd_span=$[$mjd_e-$mjd_s]

    if [ $mjd_span -lt 0 ]; then
        echo -e "$MSGERR illegal time span: from $mjd_s to $mjd_e"
        exit 1
    elif [ $mjd_span -eq   0 ]; then
        local work_dir="$proj_dir/$ymd_s/$doy_s"
    elif [ $mjd_span -lt 108 ]; then
        local work_dir="$proj_dir/$ymd_s/$doy_s-$doy_e"
    elif [ $mjd_span -ge 108 ]; then
        echo -e "$MSGERR time span too long (> 108 days): from $mjd_s to $mjd_e"
        exit 1
    fi

    mkdir -p "$work_dir" && cd "$work_dir"
    if [ $? -eq 0 ]; then
        ProcessSingleSession "$rnxo_path" "$ctrl_file" "$date_s" "$hour_s" "$date_e" "$hour_e" "$AR" \
            || echo -e "$MSGERR from $ymd_s $doy_s to $ymd_e $doy_e processing failed"
    else
        echo -e "$MSGERR no such directory: $work_dir"
    fi
}

ParseCmdArgs() { # purpose : parse command line into arguments
                 # usage   : ParseCmdArgs "$@"
    if [ $# -le 0 ]; then
        PRIDE_PPPAR_HELP
        >&2 echo ""
        PRIDE_PPPAR_INFO
        exit 1
    fi

    readonly local SITE_REGEX="^[[:alpha:]0-9]{4}$"
    readonly local PNUM_REGEX="^[+.]?[0-9]+([.][0-9]+)?$"

    local i s t iarg carg time_sec avail_num sys_num sys
    local rnxo_path rnxo_name rinex_dir ctrl_path ctrl_file
    local ymd_s hms_s ymd_e hms_e site mode plen interval freq_cmb AR
    local avail_sys edt_opt rck_opt ztd_opt htg_opt ion_opt tide_mask lam_opt pco_opt vbs_opt
    local gnss_mask map_opt rckl rckp ztdl ztdp htgp eloff

    local last_arg=${@: -1}
    case $last_arg in
    -V | --version )
        PRIDE_PPPAR_INFO && exit 1 ;;
    -H | --help )
        PRIDE_PPPAR_HELP && exit 1 ;;
    -* )
        >&2 echo -e "$MSGERR invalid argument (the last argument should be observation file): $last_arg"
        >&2 echo -e "$MSGINF use ‘pdp3 -H’ or ‘pdp3 --help’ for more information"
        exit 1
    esac

    # Parse path of observation file
    if [ -e $last_arg ]; then
        rnxo_path=$(rdlk "$last_arg")
        rnxo_name=$(basename "$rnxo_path")
        rinex_dir=$(dirname  "$rnxo_path")
        if [[ $(head -1 $last_arg | cut -c 21-21) != "O" ]]; then
            >&2 echo -e "$MSGERR unsupported RINEX observation type: $last_arg"
            >&2 echo -e "  $(head -1 $last_arg)"
            exit 1
        fi
    else
        >&2 echo -e "$MSGERR RINEX observation file doesn't exist: $last_arg"
        exit 1
    fi

    # Parse other options
    for iarg in $(seq 1 $[$#-1]); do
        case $1 in
        -?(-)+([-[:alnum:]_]) )
            case $1 in
            ## Version & Help
            -V | --version )
                PRIDE_PPPAR_INFO && exit 1
                ;;
            -H | --help )
                PRIDE_PPPAR_HELP && exit 1
                ;;
            ## Time setting
            -s | --start )
                [ -z "$ymd_s" ] && [ -z "$hms_s" ]              || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                local time=($(echo $2 | tr '/:-' ' '))
                case ${#time[@]} in
                2 ) ymd_s=$(ydoy2ymd ${time[@]}  | awk '{printf("%04d-%02d-%02d\n",$1,$2,$3)}') ;;
                3 ) ymd_s=$(echo    "${time[@]}" | awk '{printf("%04d-%02d-%02d\n",$1,$2,$3)}') ;;
                * ) throw_invalid_arg "start date" "$2" ;;
                esac
                shift 1
                if check_optional_arg "$2" "$last_arg"; then
                    local time=($(echo $2 | tr '/:-' ' '))
                    [ ${#time[@]} -eq 3 ] \
                        && hms_s=$(echo "${time[@]}" | awk '{printf("%02d:%02d:%05.2f\n",$1,$2,$3)}') \
                        && [ ${time[0]%.*} -ge 0 -a ${time[0]%.*} -le 23 ] \
                        && [ ${time[1]%.*} -ge 0 -a ${time[1]%.*} -le 59 ] \
                        && [ ${time[2]%.*} -ge 0 -a ${time[2]%.*} -le 59 ] \
                        || throw_invalid_arg "start time" "$2"
                    shift 1
                fi
                ;;
            -e | --end )
                [ -z "$ymd_e" ] && [ -z "$hms_e" ]              || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                local time=($(echo $2 | tr '/:-' ' '))
                case ${#time[@]} in
                2 ) ymd_e=$(ydoy2ymd ${time[@]}  | awk '{printf("%04d-%02d-%02d\n",$1,$2,$3)}') ;;
                3 ) ymd_e=$(echo    "${time[@]}" | awk '{printf("%04d-%02d-%02d\n",$1,$2,$3)}') ;;
                * ) throw_invalid_arg "end date" "$2" ;;
                esac
                shift 1
                if check_optional_arg "$2" "$last_arg"; then
                    local time=($(echo $2 | tr '/:-' ' '))
                    [ ${#time[@]} -eq 3 ] \
                        && hms_e=$(echo "${time[@]}" | awk '{printf("%02d:%02d:%05.2f\n",$1,$2,$3)}') \
                        && [ ${time[0]%.*} -ge 0 -a ${time[0]%.*} -le 23 ] \
                        && [ ${time[1]%.*} -ge 0 -a ${time[1]%.*} -le 59 ] \
                        && [ ${time[2]%.*} -ge 0 -a ${time[2]%.*} -le 59 ] \
                        || throw_invalid_arg "end time" "$2"
                    shift 1
                fi
                ;;
            ## General setting
            -cfg | --config )
                [ -z "$ctrl_path" ]                             || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                if [ -e "$2" ]; then
                    ctrl_path="$(rdlk "$2")"
                else
                    >&2 echo -e "$MSGERR PRIDE PPP-AR configuration file doesn't exist: $2"
                    exit 1
                fi
                shift 1
                ;;
            -sys | --system )
                [ -z "$avail_sys" ]                             || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                carg=$(echo "$2" | tr 'a-z' 'A-Z')
                avail_sys=($(sed "s/C/23/;s/./& /g" <<< "$carg"))
                gnss_mask=("G" "R" "E" "2" "3" "J")
                for s in ${avail_sys[@]}; do
                    case ${s} in
                    @(G|R|E|2|3|J) ) gnss_mask=("${gnss_mask[@]/$s}");;
                    * ) throw_invalid_arg "GNSS" "$s" ;;
                    esac
                done
                shift 1
                ;;
            -frq | --frequency )
                [ -z "$freq_cmb" ]                              || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                while true; do
                    sys_num=($(sed "s/./& /g" <<< "$2"))        || throw_invalid_arg "frequency number" "$2"
                    [ "${#sys_num[@]}" -eq "3" ]                || throw_invalid_arg "frequency number" "$2"
                    sys=$(echo "${sys_num[0]}" | tr 'a-z' 'A-Z')
                    case "$sys" in
                    "G" ) avail_num="125"    ;;
                    "R" ) avail_num="12"     ;;
                    "E" ) avail_num="15678"  ;;
                    "C" ) avail_num="125678" ;;
                    "J" ) avail_num="1256"   ;;
                     *  ) throw_invalid_arg "frequency number" "$2" ;;
                    esac
                    echo "${freq_cmb[@]}" | grep -q "$sys"      && throw_conflict_opt "$2"
                    echo "$avail_num" | grep -q "${sys_num[1]}" || throw_invalid_arg "frequency number" "$2"
                    echo "$avail_num" | grep -q "${sys_num[2]}" || throw_invalid_arg "frequency number" "$2"
                    freq_cmb+=($(echo "${sys_num[@]}" | sed "s/ //g" | tr 'a-z' 'A-Z'))
                    shift 1
                    check_optional_arg "$2" "$last_arg" || break
                done
                freq_cmb="${freq_cmb[@]}"
                ;;
            -n | --site )
                [ -z "$site" ]                                  || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                [[ "$2" =~ $SITE_REGEX ]] && site="$2"          || throw_invalid_arg "site name" "$2"
                shift 1
                ;;
            -m | --mode )
                [ -z "$mode" ]                                  || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                carg=$(echo "$2" | cut -c 1 | tr 'a-z' 'A-Z')   || throw_invalid_arg "positioning mode" "$2"
                case ${carg} in
                @(S|K|F|L) )
                    mode="$carg" ;;
                "P" )
                    mode="$carg"
                    plen=${2:1}
                    [ -n "$plen" ] || plen="300"
                    if [[ $plen =~ $PNUM_REGEX ]]                && \
                       [[ $(echo "1 <= $plen" | bc) -eq 1 ]]; then
                        plen=$(printf "%d" $[10#$plen])
                        mode="$mode:$plen"
                    else
                        throw_invalid_arg "position piece length" "$plen"
                    fi
                    ;;
                * ) throw_invalid_arg "positioning mode" "$2" ;;
                esac
                shift 1
                ;;
            -i | --interval )
                [ -z "$interval" ]                              || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                if [[ $2 =~ $PNUM_REGEX ]]                && \
                   [[ $(echo "0.02  <= $2" | bc) -eq 1 ]] && \
                   [[ $(echo "$2 <= 300.0" | bc) -eq 1 ]]; then
                    interval="$2"
                else
                    throw_invalid_arg "interval" "$2"
                fi
                shift 1
                ;;
            -f | --float )
                [ -z "$AR" ] && AR="N"                          || throw_conflict_opt "$1"
                ;;
            ## Advanced settings
            -aoff | --wapc-off )
                [ -z "$pco_opt" ]                               || throw_conflict_opt "$1"
                pco_opt="NO"
                ;;
            -c | --cutoff-elev )
                [ -z "$eloff" ]                                 || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                if [[ $2 =~ $PNUM_REGEX ]]               && \
                   [[ $(echo "0.00 <= $2" | bc) -eq 1 ]] && \
                   [[ $(echo "$2 <= 60.0" | bc) -eq 1 ]]; then
                    eloff="$2"
                else
                    throw_invalid_arg "cutoff elevation" "$2"
                fi
                shift 1
                ;;
            -l | --loose-edit )
                [ -z "$edt_opt" ]                               || throw_conflict_opt "$1"
                edt_opt="NO"
                ;;
            -hion | --high-ion )
                [ -z "$ion_opt" ]                               || throw_conflict_opt "$1"
                ion_opt="YES"
                ;;
            -h | --htg )
                [ -z "$htg_opt" ]                               || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                carg=$(echo "$2" | cut -c 1 | tr 'a-z' 'A-Z')   || throw_invalid_arg "HTG model" "$2"
                case ${carg} in
                "P" )
                    htg_opt="PWC" && htgl=${2:1} && htgl=${htgl##*:}
                    [ -n "$htgl" ] || htgl="720"
                    if [[ $htgl =~ $PNUM_REGEX ]]              && \
                       [[ $(echo "60 <= $htgl" | bc) -eq 1 ]]; then
                        htg_opt="$htg_opt:$htgl"
                    else
                        throw_invalid_arg "HTG piece length" "$htgl"
                    fi
                    ;;
                "S" )
                    htg_opt="STO"
                    ;;
                "N" )
                    htg_opt="NON"
                    ;;
                 *  )
                    throw_invalid_arg "HTG model" "$2"
                esac
                shift 1
                if check_optional_arg "$2" "$last_arg"; then
                    if [[ $2 =~ $PNUM_REGEX ]]               && \
                       [[ $(echo "0.00 <= $2" | bc) -eq 1 ]] && \
                       [[ $(echo "$2 <= 10.0" | bc) -eq 1 ]]; then
                        htgp="$2"
                    else
                        throw_invalid_arg "HTG process noise" "$2"
                    fi
                    shift 1
                fi
                ;;
            -p | --mapping-func )
                [ -z "$map_opt" ]                               || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                map_opt=$(echo $2 | tr 'a-z' 'A-Z')
                case $map_opt in
                "G" |         "GMF" ) map_opt="GMF" ;;
                "N" | "NIE" | "NMF" ) map_opt="NIE" ;;
                "1" | "V1"  | "VM1" ) map_opt="VM1" ;;
                "3" | "V3"  | "VM3" ) map_opt="VM3" ;;
                * ) throw_invalid_arg "mapping function" "$2"
                esac
                shift 1
                ;;
            -r | --rck )
                [ -z "$rck_opt" ]                               || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                carg=$(echo "$2" | cut -c 1 | tr 'a-z' 'A-Z')   || throw_invalid_arg "clock model" "$2"
                case $carg in
                "S" )
                    rck_opt="STO"
                    ;;
                "W" )
                    rck_opt="WNO"
                    ;;
                 *  )
                    throw_invalid_arg "clock model" "$2"
                esac
                shift 1
                if check_optional_arg "$2" "$last_arg"; then
                    if [[ $2 =~ $PNUM_REGEX ]]               && \
                       [[ $(echo "0.00 <= $2" | bc) -eq 1 ]] && \
                       [[ $(echo "$2 <= 10.0" | bc) -eq 1 ]]; then
                        rckp="$2"
                    else
                        throw_invalid_arg "clock process noise" "$2"
                    fi
                    shift 1
                fi
                ;;
            -toff | --tide-off )
                [ -z "$tide_mask" ]                             || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                carg=$(echo "$2" | tr 'a-z' 'A-Z')
                tide_mask=($(sed "s/./& /g" <<< "$carg"))
                for t in ${tide_mask[@]}; do
                    case ${t} in
                    @(S|O|P) ) continue ;;
                    * ) throw_invalid_arg "tide model" "$t" ;;
                    esac
                done
                shift 1
                ;;
            -v | --verbose )
                vbs_opt="YES"
                ;;
            -x | --fix-method )
                [ -z "$lam_opt" ]                               || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                case ${2} in
                "1" ) lam_opt="NO"  ;;
                "2" ) lam_opt="YES" ;;
                * ) throw_invalid_arg "fixing method" "$2"
                esac
                shift 1
                ;;
            -z | --ztd )
                [ -z "$ztd_opt" ]                               || throw_conflict_opt "$1"
                check_optional_arg "$2" "$last_arg"             || throw_require_arg  "$1"
                carg=$(echo "$2" | cut -c 1 | tr 'a-z' 'A-Z')   || throw_invalid_arg "ZTD model" "$2"
                case ${carg} in
                "P" )
                    ztd_opt="PWC" && ztdl=${2:1} && ztdl=${ztdl##*:}
                    [ -n "$ztdl" ] || ztdl="60"
                    if [[ $ztdl =~ $PNUM_REGEX ]]              && \
                       [[ $(echo " 1 <= $ztdl" | bc) -eq 1 ]]; then
                        ztd_opt="$ztd_opt:$ztdl"
                    else
                        throw_invalid_arg "ZTD piece length" "$ztdl"
                    fi
                    ;;
                "S" )
                    ztd_opt="STO"
                    ;;
                "N" )
                    ztd_opt="NON"
                    ;;
                 *  )
                    throw_invalid_arg "ZTD model" "$2"
                esac
                shift 1
                if check_optional_arg "$2" "$last_arg"; then
                    if [[ $2 =~ $PNUM_REGEX ]]               && \
                       [[ $(echo "0.00 <= $2" | bc) -eq 1 ]] && \
                       [[ $(echo "$2 <= 10.0" | bc) -eq 1 ]]; then
                        ztdp="$2"
                    else
                        throw_invalid_arg "ZTD process noise" "$2"
                    fi
                    shift 1
                fi
                ;;
            ## End
            * )
                [[ $1 == $last_arg ]] && break
                throw_invalid_opt "$1"
                ;;
            esac
            shift 1
            ;;
        "" )
            break
            ;;
        ** )
            [[ $1 == $last_arg ]] && break
            throw_invalid_opt "$1"
            ;;
        esac
    done

    # Use default config file
    local config_template_path="$(dirname $(which pdp3))/config_template"
    [ -n "$ctrl_path" ] || ctrl_path="$config_template_path"

    if [ ! -e "$ctrl_path" ]; then
        >&2 echo -e "$MSGERR PRIDE PPP-AR configuration file doesn't exist: $ctrl_path"
        exit 1
    fi

    local opt_lin=$(grep "^ .... [A-Z] " "$ctrl_path")
    local opt_num=$(echo "$opt_lin" | wc -l)
    if [ $opt_num -ne 1 ]; then
        [ $opt_num -eq 0 ] && >&2 echo -e "$MSGERR no option line to be processed: $ctrl_path"
        [ $opt_num -gt 1 ] && >&2 echo -e "$MSGERR more than one option line to be processed: $ctrl_path"
        exit 1
    fi

    # Create temporary config file
    ctrl_file=$(mktemp -u | sed "s/tmp\./config\./")
    cp -f "$ctrl_path" "$ctrl_file" && chmod 644 "$ctrl_file"
    if [ $? -ne 0 ]; then
        >&2 echo -e "$MSGERR failed to create temporary config file: $ctrl_file"
        exit 1
    fi

    # Try getting frequency combinations from the config file
    if [ -z "$freq_cmb" ]; then
        freq_cmb=$(get_ctrl "$ctrl_file" "Frequency combination")
        echo "$freq_cmb" | grep -q "WARNING" && freq_cmb=""
        [ -n "$freq_cmb" ] && [ "$freq_cmb" != "Default" ] || freq_cmb="G12 R12 E15 C26 J12"
    fi

    freq_cmb=($(echo "$freq_cmb"))
    for carg in ${freq_cmb[@]}; do
        sys_num=($(sed "s/./& /g" <<< "$carg"))
        if [ "${#sys_num[@]}" -ne "3" ]; then
            >&2 echo -e "$MSGERR invalid frequency combination: $carg"
            exit 1
        fi
        sys=$(echo "${sys_num[0]}" | tr 'a-z' 'A-Z')
        case "$sys" in
        "G" ) avail_num="125"    ;;
        "R" ) avail_num="12"     ;;
        "E" ) avail_num="15678"  ;;
        "C" ) avail_num="125678" ;;
        "J" ) avail_num="1256"   ;;
         *  )
            >&2 echo -e "$MSGERR invalid GNSS for specifying frequency combination: $carg"
            exit 1
            ;;
        esac
        for i in $(seq 1 2); do
            if ! echo "$avail_num" | grep -q "${sys_num[$i]}"; then
                >&2 echo -e "$MSGERR invalid frequency combination: $carg"
                exit 1
            fi
        done
    done

    freq_cmb="${freq_cmb[@]}"
    for sys in "G" "R" "E" "C" "J"; do
        num=$(echo "$freq_cmb" | grep -o "$sys" | wc -l)
        if [ "$num" -ge 2 ]; then
            >&2 echo -e "$MSGERR repeated frequency combination for GNSS ($sys): $freq_cmb"
            exit 1
        elif [ "$num" -eq 0 ]; then
            case "$sys" in
            "G" ) carg="G12" ;;
            "R" ) carg="R12" ;;
            "E" ) carg="E15" ;;
            "C" ) carg="C26" ;;
            "J" ) carg="J12" ;;
            esac
            freq_cmb="$freq_cmb $carg"
        fi
    done

    if grep -q "^Frequency combination" "$ctrl_file"; then
        sedi "/^Frequency combination/s/ = .*/ = $freq_cmb/" "$ctrl_file"
    else
        sedi "/^Interval/i\\Frequency combination  = $freq_cmb" "$ctrl_file"
    fi

    # Positioning mode 
    if [ -z "$mode" ]; then
        mode=$(echo "$opt_lin" | cut -c 7-7)
        [ "$mode" == "X" ] && mode="K"
    fi

    [ -n "$mode" ] && mode=$(echo $mode | tr 'a-z' 'A-Z') || mode="K"

    # Get LEO satellie name form observation file and set site name
    if [ "$mode" == "L" -a -z "$site" ]; then
        site=$(grep "MARKER NAME" "$rnxo_path" | awk -v sep='-' '{print $1sep$2}')
        [ -n "$site" ] && site="${site:0:3}${site: -1}"
    fi

    # Default as MARKER NAME or the name of observation file
    if [ -z "$site" ]; then
        site=$(grep "MARKER NAME" "$rnxo_path" | awk '{print substr($0,0,4)}')
        if [[ ! "$site" =~ $SITE_REGEX ]]; then
            if [[ $rnxo_name =~ ^[[:alpha:]0-9]{9}_.+O\.(rnx|RNX)$ ]] || \
               [[ $rnxo_name =~ ^[[:alpha:]0-9]{4}[0-9]{3}.+\.[0-9]{2}(o|O)$ ]]; then
                site=${rnxo_name:0:4}
            fi
        fi
        if [ -z "$site" ]; then
            site=$(echo "$opt_lin" | cut -c 2-5)
            [ "$site" == "xxxx" ] && site=""
        fi
        if [[ ! "$site" =~ $SITE_REGEX ]]; then
            >&2 echo -e "$MSGERR site name not found in command-line or observation file"
            >&2 echo -e "$MSGINF please use the option ‘-n’ or ‘--site’ to input the site name"
            exit 1
        fi
    fi

    [ -n "$site" ] && site=$(echo $site | tr 'A-Z' 'a-z') || site="xxxx"

    # Default as the first epoch in the first observation file
    if [ -z "$ymd_s" ] || [ -z "$hms_s" ]; then
        time_sec=$(grep -E "^(> [ 0-9]{4} [ 0-1][0-9] | [ 0-9][0-9] [ 0-1][0-9] )" "$rnxo_path")
        local time=$(echo "$time_sec" | head -1)
        if [ -n "$time" ]; then
            if [[ $time =~ ^\> ]]; then
                [ -z "$ymd_s" ] && ymd_s=$(echo "$time" | awk '{printf("%04d-%02d-%02d\n",$2,$3,$4)}')
                [ -z "$hms_s" ] && hms_s=$(echo "$time" | awk '{printf("%02d:%02d:%010.7f\n",$5,$6,$7)}')
            else
                [ -z "$ymd_s" ] && ymd_s=$(echo "$time" | awk '{yr=$1+2000;if($1>80)yr-=100;printf("%04d-%02d-%02d\n",yr,$2,$3)}')
                [ -z "$hms_s" ] && hms_s=$(echo "$time" | awk '{printf("%02d:%02d:%010.7f\n",$4,$5,$6)}')
            fi
        else
            >&2 echo -e "$MSGERR start time not found in command-line or observation file"
            >&2 echo -e "$MSGINF please use the option ‘-s’ or ‘--start’ to input the start time"
            exit 1
        fi
    fi

    if [ -z "$ymd_e" ]; then
        local tmpfobs="$rnxo_name"
    else
        # Check RINEX OBS files for multi-day processing
        local mjd_s=$(ymd2mjd $(echo "${ymd_s[*]}" | tr '-' ' '))
        local mjd_e=$(ymd2mjd $(echo "${ymd_e[*]}" | tr '-' ' '))
        if [ "$OS" == "Darwin" ]; then
            local doy_s=$(date -j -f "%Y-%m-%d" "$date_s" +"%j")
            local doy_e=$(date -j -f "%Y-%m-%d" "$date_e" +"%j")
        else
            local doy_s=$(date -d "${ymd_s[*]}" +"%j")
            local doy_e=$(date -d "${ymd_e[*]}" +"%j")
        fi

        readonly local RNXO2D_GLOB="${rnxo_name:0:4}${doy_s}0.${ymd_s:2:2}@(o|O)"
        readonly local RNXO2H_GLOB="${rnxo_name:0:4}${doy_s}[a-x].${ymd_s:2:2}@(o|O)"
        readonly local RNXO3D_GLOB="${rnxo_name:0:9}_?_${ymd_s:0:4}${doy_s}0000_01D_???_?O.@(rnx|RNX)"
        readonly local RNXO3H_GLOB="${rnxo_name:0:9}_?_${ymd_s:0:4}${doy_s}[0-9][0-9]00_01H_???_?O.@(rnx|RNX)"

        local mjd hh hh_s hh_e al
        local tmpfobs tmpydoy
        for mjd in $(seq $mjd_s $mjd_e); do
            tmpydoy=($(mjd2ydoy $mjd))
            case "$rnxo_name" in
            $RNXO2D_GLOB )
                tmpfobs="${rnxo_name:0:4}${tmpydoy[1]}0.${tmpydoy[0]:2:2}${rnxo_name:11}"
                [ -f "$rinex_dir/$tmpfobs" ] || >&2 echo -e "$MSGWAR $rinex_dir/$tmpfobs doesn't exist"
                ;;
            $RNXO2H_GLOB )
                [ "$mjd" -eq "$mjd_s" ] && [ -n "$hms_s" ] && hh_s="${hms_s:0:2}" || hh_s="00"
                [ "$mjd" -eq "$mjd_e" ] && [ -n "$hms_e" ] && hh_e="${hms_e:0:2}" || hh_e="23"
                for hh in $(seq -w $hh_s $hh_e); do
                    al=$(printf "\\x$(printf "%x" $[$(printf "%d" "'a")+$[10#$hh]])")
                    tmpfobs="${rnxo_name:0:4}${tmpydoy[1]}${al}.${tmpydoy[0]:2:2}${rnxo_name:11}"
                    [ -f "$rinex_dir/$tmpfobs" ] || >&2 echo -e "$MSGWAR $rinex_dir/$tmpfobs doesn't exist"
                done
                ;;
            $RNXO3D_GLOB )
                tmpfobs="${rnxo_name:0:12}${tmpydoy[0]}${tmpydoy[1]}${rnxo_name:19}"
                [ -f "$rinex_dir/$tmpfobs" ] || >&2 echo -e "$MSGWAR $rinex_dir/$tmpfobs doesn't exist"
                ;;
            $RNXO3H_GLOB )
                [ "$mjd" -eq "$mjd_s" ] && [ -n "$hms_s" ] && hh_s="${hms_s:0:2}" || hh_s="00"
                [ "$mjd" -eq "$mjd_e" ] && [ -n "$hms_e" ] && hh_e="${hms_e:0:2}" || hh_e="23"
                for hh in $(seq -w $hh_s $hh_e); do
                    tmpfobs="${rnxo_name:0:12}${tmpydoy[0]}${tmpydoy[1]}${hh}${rnxo_name:21}"
                    [ -f "$rinex_dir/$tmpfobs" ] || >&2 echo -e "$MSGWAR $rinex_dir/$tmpfobs doesn't exist"
                done
                ;;
            * )
                tmpfobs="$rnxo_name"
                >&2 echo -e "$MSGWAR illegal naming convention: $tmpfobs"
                >&2 echo -e "$MSGINF please make ensure the RINEX observation file contains enough observation data to be processed"
                break
                ;;
            esac
        done
    fi

    # Default as the last epoch in the last observation file
    if [ -z "$ymd_e" ] || [ -z "$hms_e" ]; then
        time_sec=$(grep -E "^(> [0-9]{4} [ 0-1][0-9] | [ 0-9][0-9] [ 0-1][0-9] )" "$rinex_dir/$tmpfobs")
        local time=$(echo "$time_sec" | tail -1)
        if [ -n "$time" ]; then
            if [[ $time =~ ^\> ]]; then
                [ -z "$ymd_e" ] && ymd_e=$(echo "$time" | awk '{printf("%04d-%02d-%02d\n",$2,$3,$4)}')
                [ -z "$hms_e" ] && hms_e=$(echo "$time" | awk '{printf("%02d:%02d:%010.7f\n",$5,$6,$7)}')
            else
                [ -z "$ymd_e" ] && ymd_e=$(echo "$time" | awk '{yr=$1+2000;if($1>80)yr-=100;printf("%04d-%02d-%02d\n",yr,$2,$3)}')
                [ -z "$hms_e" ] && hms_e=$(echo "$time" | awk '{printf("%02d:%02d:%010.7f\n",$4,$5,$6)}')
            fi
        else
            >&2 echo -e "$MSGERR end time not found in command-line or observation file"
            >&2 echo -e "$MSGINF please use the option ‘-e’ or ‘--end’ to input the end time"
            exit 1
        fi
    fi

    # Check time span
    if [ "$OS" == "Darwin" ]; then
        local sec_s=$(date -j -f "%Y-%m-%d %H:%M:%S" "${ymd_s} ${hms_s%.*}" +"%s" | awk '{printf("%.3f", $1 + ("0." "'${hms_s#*.}'"))}')
        local sec_e=$(date -j -f "%Y-%m-%d %H:%M:%S" "${ymd_e} ${hms_e%.*}" +"%s" | awk '{printf("%.3f", $1 + ("0." "'${hms_s#*.}'"))}')
    else
        local sec_s=$(date -d "$ymd_s $hms_s" +"%s.%3N")
        local sec_e=$(date -d "$ymd_e $hms_e" +"%s.%3N")
    fi
    local sspan=$(echo "$sec_e - $sec_s" | bc)
    if [[ $(echo "$sspan <= 0" | bc) -eq 1 ]]; then
        >&2 echo -e "$MSGERR illegal time span: from $ymd_s $hms_s to $ymd_e $hms_e"
        exit 1
    fi

    local session_time="${ymd_s[@]//-/ } ${hms_s[@]//:/ } ${sspan}"
    sedi "/^Session time/s/ = .*/ = $session_time/" "$ctrl_file"

    # Try getting observation interval option from the config file
    [ -n "$interval" ] || interval=$(get_ctrl "$ctrl_file" "Interval")
    [ -z "$time_sec" ] && time_sec=$(grep -E "^(> [0-9]{4} [ 0-1][0-9] | [ 0-9][0-9] [ 0-1][0-9] )" "$rnxo_path")
    local obsintvl=$(echo "$time_sec" | awk 'BEGIN {
                                                 mdif = 30
                                             } {
                                                 if ($1 == ">") {
                                                     this_sec = $5*3600+$6*60+$7
                                                 } else {
                                                     this_sec = $4*3600+$5*60+$6
                                                 }
                                                 if (last_sec != "") {
                                                     vdif = this_sec - last_sec
                                                     if (vdif < 0) vdif *= -1
                                                     if (vdif < mdif && vdif != 0) mdif = vdif
                                                 }
                                                 last_sec = this_sec
                                             } END {
                                                 print(mdif)
                                             }')

    if [[ -n "$interval" ]] && [[ "$interval" != "Default" ]]; then
        if [[ $(echo "$interval < $obsintvl" | bc) -eq 1 ]]; then
            >&2 echo -e "$MSGERR input interval is shorter than observation interval: $interval < $obsintvl"
            exit 1
        fi
    else
        ## Align to the nearest candidate
        local last_can last_dif this_dif
        local cand=("86400" "300" "60" "30" "25" "20" "15" "10" "5" "2" "1" "0.5" "0.25" "0.2" "0.1" "0.05" "0.02" "-86400")
        for i in $(seq 1 $[${#cand[@]}-1]); do
            last_can=${cand[$[$i-1]]}
            last_dif=$(echo "$obsintvl" | awk '{print("'${last_can}'"-$0)}')
            this_dif=$(echo "$obsintvl" | awk '{print($0-"'${cand[$i]}'")}')
            if [[ $(echo "$last_dif < $this_dif" | bc) -eq 1 ]]; then
                if [[ $(echo "$obsintvl == $last_can" | bc) -ne 1 ]]; then
                    >&2 echo -e "$MSGWAR observation interval rounded to the nearest candidate: $obsintvl -> $last_can"
                fi
                obsintvl="$last_can" && break
            fi
        done
        interval="$obsintvl"
    fi

    if [[ $(echo "0.02 > $interval" | bc) -eq 1 ]]; then
       >&2 echo -e "$MSGWAR observation interval too small, rounded to the minimum candidate: $interval -> 0.02"
       interval="0.02"
    fi

    if [[ $(echo "$interval > 300.0" | bc) -eq 1 ]]; then
       >&2 echo -e "$MSGWAR observation interval too large, rounded to the maximum candidate: $interval -> 300.0"
       interval="300.0"
    fi

    sedi "/^Interval/s/ = .*/ = $interval/" "$ctrl_file"

    # Piece length for PWC mode
    if [ -z "$plen" ]; then
        plen=$(echo "$opt_lin" | awk '{print($13)}')
        [[ ! $plen =~ $PNUM_REGEX ]] && plen="300"
    fi

    if [ "${mode:0:1}" == "P" ]; then
        if [[ $(echo "$plen < $interval" | bc) -eq 1 ]]; then
            >&2 echo -e "$MSGERR piece-length shorter than observation interval: $plen < $interval"
            exit 1
        fi
        if [[ $(echo "$plen > $sspan" | bc) -eq 1 ]]; then
            >&2 echo -e "$MSGWAR piece-length longer than observation period: $plen > $sspan"
        fi
    fi

    # Editing mode
    [ -n "$edt_opt" ] || edt_opt=$(get_ctrl "$ctrl_file" "Strict editing")
    [ "$edt_opt" == "Default" ] && edt_opt="YES"
    sedi "/^Strict editing/s/ = .*/ = $edt_opt/" "$ctrl_file"

    [ "$edt_opt" == "YES" ] && min_sspan="600.0" || min_sspan="30.0"
    if [[ $(echo "$sspan < $min_sspan" | bc) -eq 1 ]]; then
        >&2 echo -e "$MSGERR observation period too short: $sspan < $min_sspan"
        exit 1
    fi

    # RCK model
    [ -n "$rck_opt" ] || rck_opt=$(get_ctrl "$ctrl_file" "RCK model")
    [ "$rck_opt" == "Default" ] && rck_opt="WNO"

    sedi "/^RCK model/s/ = .*/ = $rck_opt/" "$ctrl_file"

    if [ -z "$rckp" ]; then
        local rckp=$(echo "$opt_lin" | awk '{print($5)}')
        if [[ ! $rckp =~ $PNUM_REGEX ]]; then
            case ${rck_opt:0:3} in
            "STO" ) rckp="0.001" ;;
            "WNO" ) rckp="0.000" ;;
              *   ) rckp="0.000" ;;
            esac
        fi
    fi

    # ZTD model
    [ -n "$ztd_opt" ] || ztd_opt=$(get_ctrl "$ctrl_file" "ZTD model")

    if [ "$mode" == "L" -a "$ztd_opt" != "NON" ]; then
        [ "$ztd_opt" != "Default" ] && \
            >&2 echo -e "$MSGWAR disable ZTD model for LEO satellite: $ztd_opt -> NON"
        ztd_opt="NON"
    fi

    [ "$ztd_opt" == "Default" ] && ztd_opt="STO"

    sedi "/^ZTD model/s/ = .*/ = $ztd_opt/" "$ctrl_file"

    if [ -z "$ztdp" ]; then
        local ztdp=$(echo "$opt_lin" | awk '{print($8)}')
        if [[ ! $ztdp =~ $PNUM_REGEX ]]; then
            case ${ztd_opt:0:3} in
            "STO" ) ztdp=".0004" ;;
            "PWC" ) ztdp="0.020" ;;
              *   ) ztdp="0.020" ;;
            esac
        fi
    fi

    ztdp=$(printf "%5s" $ztdp)

    # HTG model
    [ -n "$htg_opt" ] || htg_opt=$(get_ctrl "$ctrl_file" "HTG model")

    if [ "$mode" == "L" -a "$htg_opt" != "NON" ]; then
        [ "$htg_opt" != "Default" ] && \
            >&2 echo -e "$MSGWAR disable HTG model for LEO satellite: $htg_opt -> NON"
        ztd_opt="NON"
    fi

    if [ "$htg_opt" == "Default" ]; then
        [ "$mode" == "F" -o $mode == "S" ] && htg_opt="PWC:720" || htg_opt="NON"
    fi

    sedi "/^HTG model/s/ = .*/ = $htg_opt/" "$ctrl_file"

    if [ -z "$htgp" ]; then
        local htgp=$(echo "$opt_lin" | awk '{print($10)}')
        if [[ ! $htgp =~ $PNUM_REGEX ]]; then
            case ${htg_opt:0:3} in
            "STO" ) htgp=".0004" ;;
            "PWC" ) htgp="0.002" ;;
              *   ) htgp="0.002" ;;
            esac
        fi
    fi

    htgp=$(printf "%5s" $htgp)

    # High-order ionospheric delay model
    [ -n "$ion_opt" ] || ion_opt=$(get_ctrl "$ctrl_file" "Iono 2nd")
    [ "$ion_opt" == "Default" ] && ion_opt="NO"
    sedi "/^Iono 2nd/s/ = .*/ = $ion_opt/" "$ctrl_file"

    # Tide correction model
    local tide_mode
    if [ -z "$tide_mask" ]; then
        tide_mode=$(get_ctrl "$ctrl_file" "Tides")
        if [ "$tide_mode" == "Default" ]; then
            [ "$mode" == "L" ] && tide_mode="NON" || tide_mode="SOLID/OCEAN/POLE"
        fi
    else
        tide_mode=("SOLID" "OCEAN" "POLE")
        for t in ${tide_mask[@]}; do
            t=$(echo $t | tr 'a-z' 'A-Z')
            [ "$t" == "S" ] && tide_mode=("${tide_mode[@]/SOLID}")
            [ "$t" == "O" ] && tide_mode=("${tide_mode[@]/OCEAN}")
            [ "$t" == "P" ] && tide_mode=("${tide_mode[@]/POLE}")
        done
        tide_mode=$(echo "${tide_mode[@]}" | sed "s/^ *//; s/ *$//; s/  */\//g")
        [ -n "$tide_mode" ] || tide_mode="NON"
    fi

    if [ "$mode" == "L" -a "$tide_mode" != "NON" ]; then
        >&2 echo -e "$MSGWAR disable tidal correction for LEO satellite: $tide_mode -> NON"
        tide_mode="NON"
    fi

    sedi "/^Tides/s/ = .*/ = ${tide_mode//\//\\/}/" "$ctrl_file"

    # Ambiguity resolution
    [ -n "$AR" ] || AR="A"

    [ -n "$lam_opt" ] || lam_opt=$(get_ctrl "$ctrl_file" "Ambiguity co-var")
    if [ "$lam_opt" == "Default" ]; then
        ## max observation time set for LAMBDA is 6 hours
        [[ $(echo "$sspan <= 21600.0" | bc) -eq 1 ]] && lam_opt="YES"     || lam_opt="NO"
        [ "$mode" == "L" ] && [[ $(echo "$sspan >= 3600.0" | bc) -eq 1 ]] && lam_opt="NO"
    fi

    sedi "/^Ambiguity co-var/s/ = .*/ = $lam_opt/" "$ctrl_file"

    # PCO on wide-lane (default as YES)
    [ -n "$pco_opt" ] || pco_opt=$(get_ctrl "$ctrl_file" "PCO on wide-lane" | tr 'a-z' 'A-Z')
    [[ "$pco_opt" == "NO" ]] || pco_opt="YES"
    sedi "/^PCO on wide-lane/s/ = .*/ = $pco_opt/" "$ctrl_file"

    # Verbose output (default as NO)
    [ -n "$vbs_opt" ] || vbs_opt=$(get_ctrl "$ctrl_file" "Verbose output" | tr 'a-z' 'A-Z')
    [[ "$vbs_opt" == "YES" ]] || vbs_opt="NO"
    sedi "/^Verbose output/s/ = .*/ = $vbs_opt/" "$ctrl_file"

    # GNSS
    for s in ${gnss_mask[@]}; do
        s=$(echo $s | tr 'a-z' 'A-Z')
        case $s in
        "2" ) prn_mask=($(seq -f  "C%02g"  1 16)) ;;
        "3" ) prn_mask=($(seq -f  "C%02g" 17 99)) ;;
         *  ) prn_mask=($(seq -f "$s%02g"  1 99)) ;;
        esac
        for prn in ${prn_mask[@]}; do
            sedi "/^ $prn /s/^ /#/" "$ctrl_file"
        done
    done

    # Disable ambiguity resolution when process with GLONASS only
    grep -q "^ [GECJ][0-9][0-9] " "$ctrl_file" || AR="N"

    # Mapping function
    [ -n "$map_opt" ] || map_opt=$(echo "$opt_lin" | awk '{print($3)}')

    if [ "$mode" == "L" -a "$map_opt" != "NON" ]; then
        [ "$map_opt" != "Default" ] && \
            >&2 echo -e "$MSGWAR disable mapping function for LEO satellite: $map_opt -> NON"
        map_opt="NON"
    fi

    [ "$map_opt" == "XXX" ] && map_opt="GMF"

    # Cutoff elevation
    [ -n "$eloff" ] || eloff=$(echo "$opt_lin" | awk '{print($6)}')
    if [[ ! $eloff =~ $PNUM_REGEX ]]; then
        [ "$mode" == "L" ] && eloff="0" || eloff="7"
    fi

    eloff=$(echo "$eloff" | awk '{printf("%2d",$0)}')

    # Modify option line
    local clkm=$(echo "$opt_lin" | awk '{print($4)}')
    local ztdm=$(echo "$opt_lin" | awk '{print($7)}')
    local htgm=$(echo "$opt_lin" | awk '{print($9)}')
    local ragm=$(echo "$opt_lin" | awk '{print($11)}')
    local phsc=$(echo "$opt_lin" | awk '{print($12)}')
    local poxm=$(echo "$opt_lin" | awk '{print($14)}')
    local poym=$(echo "$opt_lin" | awk '{print($15)}')
    local pozm=$(echo "$opt_lin" | awk '{print($16)}')

    opt_lin=" $site ${mode:0:1}  $map_opt $clkm $rckp $eloff $ztdm $ztdp $htgm $htgp $ragm $phsc $plen $poxm $poym $pozm"
    sedi "s/^ .... [A-Z] .*/$opt_lin/" "$ctrl_file"

    # Return
    echo "$rnxo_path"
    echo "$ctrl_path"
    echo "$ctrl_file"
    echo "$ymd_s"
    echo "${hms_s:0:12}"
    echo "$ymd_e"
    echo "${hms_e:0:12}"
    echo "$mode"
    echo "$AR"
}

check_optional_arg() { # purpose : check if optional argument is existing
                       # usage   : check_optional_arg this_arg last_arg
    local this_arg=$1
    local last_arg=$2
    [[ -z $this_arg ]] && return 1
    [[ $this_arg =~ ^-{1,2}   ]] && return 1
    [[ $this_arg == $last_arg ]] && return 1
    return 0
}

throw_conflict_opt() { # purpose : throw exception message and exit when option conflicts with a previous option
                       # usage   : throw_invalid_arg opt
    local opt=$1
    >&2 echo "$SCRIPT_NAME: conflicting option '$opt'"
    >&2 echo "Try '$SCRIPT_NAME --help' for more information."
    exit 1
}

throw_invalid_arg() { # purpose : throw exception message and exit when option got an invalid argument
                      # usage   : throw_invalid_arg optlable argument
    local optlable=$1
    local argument=$2
    >&2 echo "$SCRIPT_NAME: invalid $optlable: ‘$argument’"
    >&2 echo "Try '$SCRIPT_NAME --help' for more information."
    exit 1
}

throw_invalid_opt() { # purpose : throw exception message and exit when an invalid option occurs
                      # usage   : throw_invalid_opt opt
    local opt=$1
    case $opt in
    --+([-[:alnum:]_]) ) local detail="unrecognized option '${opt}'" ;;
     -+([-[:alnum:]_]) ) local detail="invalid option -- '${opt:1}'" ;;
      * )                local detail="invalid argument -- '${opt}'" ;;
    esac
    >&2 echo "$SCRIPT_NAME: $detail"
    >&2 echo "Try '$SCRIPT_NAME --help' for more information."
    exit 1
}

throw_require_arg() { # purpose : throw exception message and exit when option did not get its argument
                      # usage   : throw_require_arg opt
    local opt=$1
    case $opt in
    --+([-[:alnum:]_]) ) local detail="option '${opt:0}' requires an argument"    ;;
     -+([-[:alnum:]_]) ) local detail="option requires an argument -- '${opt:1}'" ;;
      * )                local detail="invalid argument -- '${opt:0}'"            ;;
    esac
    >&2 echo "$SCRIPT_NAME: $detail"
    >&2 echo "Try '$SCRIPT_NAME --help' for more information."
    exit 1
}

CheckExecutables() { # purpose : check whether all needed executables are callable
                     # usage   : CheckExecutables
    echo -e "$MSGSTA CheckExecutables ..."
    for exceu in "arsig" "get_ctrl" "lsq" "redig" "sp3orb" "spp" "tedit"; do
        if ! which $exceu > /dev/null 2>&1; then
            echo -e "$MSGERR PRIDE PPP-AR executable file $exceu not found"
            return 1
        fi
    done
    for exceu in "leoatx.py" "merge2brdm.py" "pbopos" "pso2kin.py" "xyz2enu"; do
        if ! which $exceu > /dev/null 2>&1; then
            echo -e "$MSGWAR PRIDE PPP-AR executable file $exceu not found"
        fi
    done
    for exceu in "awk" "diff" "sed"; do
        if ! which $exceu > /dev/null 2>&1; then
            echo -e "$MSGERR system tool $exceu not found"
            return 1
        fi
    done
    for exceu in "curl" "gunzip" "wget"; do
        if ! which $exceu > /dev/null 2>&1; then
            echo -e "$MSGWAR system tool $exceu not found"
        fi
    done
    echo -e "$MSGSTA CheckExecutables done"
}

PRIDE_PPPAR_INFO() { # purpose : print information for PRIDE PPP-AR
                     # usage   : PRIDE_PPPAR_INFO
    >&2 echo "© GNSS Research Center of Wuhan University, 2023"
    >&2 echo "  GNSS PPP & PPP-AR data processing with PRIDE PPP-AR version $VERSION_NUM"
}

PRIDE_PPPAR_HELP() { # purpose : print usage for PRIDE PPP-AR
                     # usage   : PRIDE_PPPAR_HELP
    >&2 echo "Usage: $SCRIPT_NAME [options] <obs-file>"
    >&2 echo ""
    >&2 echo "  All character type arguments could be upper-case or lower-case"
    >&2 echo ""
    >&2 echo "Start up:"
    >&2 echo ""
    >&2 echo "  -V, --version                              display version of this script"
    >&2 echo ""
    >&2 echo "  -H, --help                                 print this help"
    >&2 echo ""
    >&2 echo "Common options:"
    >&2 echo ""
    >&2 echo "  -cfg <file>, --config <file>               configuration file for PRIDE PPP-AR 3"
    >&2 echo ""
    >&2 echo "  -sys <char>, --system <char>               GNSS to be processed, select one or more from \"GREC23J\":"
    >&2 echo "                                             -----+------------------------+-----+-------------------------"
    >&2 echo "                                               G  |  GPS                   |  R  |  GLONASS                "
    >&2 echo "                                               E  |  Galileo               |  C  |  BeiDou-2 and BeiDou-3  "
    >&2 echo "                                               2  |  BeiDou-2 only         |  3  |  BeiDou-3 only          "
    >&2 echo "                                               J  |  QZSS                  |     |                         "
    >&2 echo "                                             -----+------------------------+-----+-------------------------"
    >&2 echo "                                               * default: all GNSS"
    >&2 echo ""
    >&2 echo "  -frq <char>, --frequency <char>            frequencies to form ionosphere-free combination, select from:"
    >&2 echo "                                             -----+-------+-------+-------+-------+-------+-------+--------"
    >&2 echo "                                                  |   1   |   2   |   5   |   6   |   7   |   8   |        "
    >&2 echo "                                             -----+-------+-------+-------+-------+-------+-------+--------"
    >&2 echo "                                               G  |  L1   |  L2   |  L5   |       |       |       |        "
    >&2 echo "                                               R  |  L1   |  L2   |       |       |       |       |        "
    >&2 echo "                                               E  |  E1   |       |  E5a  |  E6   |  E5b  |  E5   |        "
    >&2 echo "                                               C  |  B1C  |  B1I  |  B2a  |  B3I  |  B2b  |  B2   |        "
    >&2 echo "                                               J  |  L1   |  L2   |  L5   |  L6   |       |       |        "
    >&2 echo "                                             -----+-------+-------+-------+-------+-------+-------+--------"
    >&2 echo "                                               input as: \"G12 R12 E15 C26 J12\" (default setting)"
    >&2 echo ""
    >&2 echo "  -m <char[length]>, --mode <char[length]>   positioning mode, select one from \"S/P/K/F/L\":"
    >&2 echo "                                             -----+--------------+-----+--------------+-----+--------------"
    >&2 echo "                                               S  |  static      |  P  | piece-wise   |  K  |  kinematic   "
    >&2 echo "                                             -----+--------------+-----+--------------+-----+--------------"
    >&2 echo "                                               F  |  fixed       |  L  | LEO sat      |     |              "
    >&2 echo "                                             -----+--------------+-----+--------------+-----+--------------"
    >&2 echo "                                               * default: kinematic mode (K)"
    >&2 echo "                                               * length argument applied for P mode only, input in seconds"
    >&2 echo "                                               * P mode default as P300 (PWC model with length as 300 s)"
    >&2 echo ""
    >&2 echo ""
    >&2 echo "  -s <date [time]>, --start <date [time]>    start date (and time) for processing, format:"
    >&2 echo "                                             --------+--------------------------+--------+-----------------"
    >&2 echo "                                               date  |  yyyy/mm/dd or yyyy/doy  |  time  |  hh:mm:ss       "
    >&2 echo "                                             --------+--------------------------+--------+-----------------"
    >&2 echo "                                               * default: first observation epoch in the obs-file"
    >&2 echo ""
    >&2 echo "  -e <date [time]>, --end <date [time]>      end date (and time) for processing, format:"
    >&2 echo "                                             --------+--------------------------+--------+-----------------"
    >&2 echo "                                               date  |  yyyy/mm/dd or yyyy/doy  |  time  |  hh:mm:ss       "
    >&2 echo "                                             --------+--------------------------+--------+-----------------"
    >&2 echo "                                               * default: last observation epoch in the obs-file"
    >&2 echo ""
    >&2 echo "  -n <char>, --site <char>                   site name for processing, format: NNNN"
    >&2 echo "                                               * default: MARKER NAME in the obs-file, or the first four"
    >&2 echo "                                                   characters of the RINEX filename"
    >&2 echo ""
    >&2 echo "  -i <num>,  --interval <num>                processing interval in seconds, 0.02 <= interval <= 300"
    >&2 echo "                                               * default: minimal observation interval in the obs-file"
    >&2 echo ""
    >&2 echo "Advanced options:"
    >&2 echo ""
    >&2 echo "  -aoff, --wapc-off                          disable APC correction on the Melbourne-Wubbena combination"
    >&2 echo ""
    >&2 echo "  -c <num>, --cutoff-elev <num>              cutoff elevation in degrees, 0 <= elevation <=60"
    >&2 echo "                                               * default: 7 degrees"
    >&2 echo ""
    >&2 echo "  -f, --float                                disable ambiguity resolution"
    >&2 echo ""
    >&2 echo "  -h <char[length] [num]>, --htg <char[length] [num]>"
    >&2 echo "                                             HTG model, piece length and process noise:"
    >&2 echo "                                                  input      model     length          process noise       "
    >&2 echo "                                             --------------+-------+-----------+---------------------------"
    >&2 echo "                                               S    .0010  |  STO  |           |      0.0010 m/sqrt(s)     "
    >&2 echo "                                               S           |  STO  |           |      0.0004 m/sqrt(s)     "
    >&2 echo "                                               P720        |  PWC  |  720 min  |      0.002  m/sqrt(h)     "
    >&2 echo "                                               P60  0.004  |  PWC  |   60 min  |      0.004  m/sqrt(h)     "
    >&2 echo "                                               P           |  PWC  |  720 min  |      0.002  m/sqrt(h)     "
    >&2 echo "                                               N           |  NON  |           |                           "
    >&2 echo "                                             --------------+-------+-----------+---------------------------"
    >&2 echo "                                               * HTG - horizontal tropospheric gradient"
    >&2 echo "                                               * STO - stochastic walk"
    >&2 echo "                                               * PWC - piece-wise constant"
    >&2 echo "                                               * NON - none"
    >&2 echo "                                               * default: PWC model with process noise as 0.002 m/sqrt(h)"
    >&2 echo "                                                   for static and fixed mode, NON for the other modes"
    >&2 echo ""
    >&2 echo "  -hion, --high-ion                          use 2nd ionospheric delay model with CODE's GIM products"
    >&2 echo ""
    >&2 echo "  -l, --loose-edit                           disable strict editing"
    >&2 echo ""
    >&2 echo "  -p <char>, --mapping-func <char>           mapping function (MF), select one from \"G/N/V1/V3\""
    >&2 echo "                                             -------+-----------------------+------+-----------------------"
    >&2 echo "                                                G   |  Global MF            |  V1  |  Vienna MF 1          "
    >&2 echo "                                                N   |  Niell MF             |  V3  |  Vienna MF 3          "
    >&2 echo "                                             -------+-----------------------+------+-----------------------"
    >&2 echo "                                               * default: global mapping function (G)"
    >&2 echo ""
    >&2 echo "  -r <char[length] [num]>, --rck <char[length] [num]>"
    >&2 echo "                                             Receiver clock model, process noise:"
    >&2 echo "                                                  input      model     length          process noise       "
    >&2 echo "                                             --------------+-------+-----------+---------------------------"
    >&2 echo "                                               S    .0040  |  STO  |           |      0.0040 m/sqrt(s)     "
    >&2 echo "                                               S           |  STO  |           |      0.0010 m/sqrt(s)     "
    >&2 echo "                                               W           |  WNO  |           |                           "
    >&2 echo "                                             --------------+-------+-----------+---------------------------"
    >&2 echo "                                               * STO - stochastic walk"
    >&2 echo "                                               * WNO - white noise"
    >&2 echo "                                               * default: WNO model (W)"
    >&2 echo ""
    >&2 echo "  -toff <char>, --tide-off <char>            disable tide correction, select one or more from \"SOP\":"
    >&2 echo "                                             -----+--------------+-----+--------------+-----+--------------"
    >&2 echo "                                               S  |  solid       |  O  | ocean        |  P  |  pole        "
    >&2 echo "                                             -----+--------------+-----+--------------+-----+--------------"
    >&2 echo "                                               * default: apply all tide corrections"
    >&2 echo ""
    >&2 echo "  -v, --verbose                              output details of ambiguity resolution"
    >&2 echo ""
    >&2 echo "  -x <num>, --fix-method <num>               ambiguity fixing method, choose 1 or 2:"
    >&2 echo "                                             -----+------------------------+-----+-------------------------"
    >&2 echo "                                               1  |  rounding              |  2  |  LAMBDA                 "
    >&2 echo "                                             -----+------------------------+-----+-------------------------"
    >&2 echo "                                               * default: rounding for long observation time (> 6 h)"
    >&2 echo "                                                          LAMBDA for short observation time (<= 6 h)"
    >&2 echo ""
    >&2 echo "  -z <char[length] [num]>, --ztd <char[length] [num]>"
    >&2 echo "                                             ZTD model, piece length and process noise:"
    >&2 echo "                                                  input      model     length          process noise       "
    >&2 echo "                                             --------------+-------+-----------+---------------------------"
    >&2 echo "                                               S    .0010  |  STO  |           |      0.0010 m/sqrt(s)     "
    >&2 echo "                                               S           |  STO  |           |      0.0004 m/sqrt(s)     "
    >&2 echo "                                               P720        |  PWC  |  720 min  |      0.02   m/sqrt(h)     "
    >&2 echo "                                               P60  0.040  |  PWC  |   60 min  |      0.04   m/sqrt(h)     "
    >&2 echo "                                               P           |  PWC  |   60 min  |      0.02   m/sqrt(h)     "
    >&2 echo "                                               N           |  NON  |           |                           "
    >&2 echo "                                             --------------+-------+-----------+---------------------------"
    >&2 echo "                                               * ZTD - zenith total delay of troposphere"
    >&2 echo "                                               * STO - stochastic walk"
    >&2 echo "                                               * PWC - piece-wise constant"
    >&2 echo "                                               * NON - none"
    >&2 echo "                                               * default: STO model with process noise as .0004 m/sqrt(s)"
    >&2 echo ""
    >&2 echo "Examples:"
    >&2 echo ""
    >&2 echo "  pdp3 abmf0010.20o                          single-day processing"
    >&2 echo ""
    >&2 echo "  pdp3 -s 2020/1 -e 2020/3 abmf0010.20o      multi-day processing from 2020/001 to 2020/003"
    >&2 echo ""
    >&2 echo "  More details refer to PRIDE PPP-AR manual and repository"
    >&2 echo "    https://github.com/PrideLab/PRIDE-PPPAR/"
}

ProcessSingleSession() { # purpose : process data of a single observation session
                         # usage   : ProcessSingleSession rnxo_path ctrl_file ymd_s hms_s ymd_e hms_e AR(A/Y/N)
    local rnxo_path="$1"
    local ctrl_file="$2"
    local ymd_s="$3"
    local hms_s="$4"
    local ymd_e="$5"
    local hms_e="$6"
    local AR="$7"

    local interval=$(get_ctrl "$ctrl_file" "Interval")
    local site=$(grep "^ .... [A-Z]" "$ctrl_file" | cut -c 2-5)
    local mode=$(grep "^ .... [A-Z]" "$ctrl_file" | cut -c 7-7)

    local rinex_dir=$(dirname  "$rnxo_path")
    local rnxo_name=$(basename "$rnxo_path")

    if [ "$OS" == "Darwin" ]; then
        local doy_s=$(date -j -f "%Y-%m-%d" "$date_s" +"%j")
        local doy_e=$(date -j -f "%Y-%m-%d" "$date_e" +"%j")
    else
        local doy_s=$(date -d "$ymd_s" +"%j")
        local doy_e=$(date -d "$ymd_e" +"%j")
    fi
    local ymd_s=($(echo "$ymd_s" | tr '-' ' '))
    local ymd_e=($(echo "$ymd_e" | tr '-' ' '))
    local mon_s=${ymd_s[1]}
    local mon_e=${ymd_e[1]}
    local day_s=${ymd_s[2]}
    local day_e=${ymd_e[2]}
    local mjd_s=$(ymd2mjd ${ymd_s[*]})
    local mjd_e=$(ymd2mjd ${ymd_e[*]})

    echo -e "$MSGSTA ProcessSingleSession from $ymd_s $doy_s to $ymd_e $doy_e ..."

    CleanAll "$ymd_s" "$doy_s"

    # Prepare config
    mv -f "$ctrl_file" . && ctrl_file=$(basename "$ctrl_file")
    if [ $? -ne 0 ]; then
        echo -e "$MSGERR temporary config file doesn't exist: $ctrl_file"
        return 1
    fi

    # Prepare tables
    local table_dir=$(get_ctrl "$ctrl_file" "Table directory" | sed "s/^[ ]*//; s/[ ]*$//; s#^~#$HOME#")
    PrepareTables "$mjd_s" "$mjd_e" "$table_dir" || return 1
    if [ $? -ne 0 ]; then
        echo -e "$MSGERR PrepareTables failed"
        return 1
    fi

    # RINEX-OBS check
    local rinexobs="$rinex_dir/$rnxo_name"
    [ ! -f "$rinexobs" ] && echo -e "$MSGWAR $rinexobs doesn't exist" && return 1

    local obs_cand sys_num sys num
    local freq_cmb=($(get_ctrl "$ctrl_file" "Frequency combination"))
    for sys_num in ${freq_cmb[@]}; do
        sys=${sys_num:0:1}
        grep -q "^ $sys[0-9][0-9] " "$ctrl_file" || continue
        obs_cand+=("${sys}XX")
        num=${sys_num:1:1}
        obs_cand+=("${sys}C${num}")
        obs_cand+=("${sys}L${num}")
        num=${sys_num:2:1}
        obs_cand+=("${sys}C${num}")
        obs_cand+=("${sys}L${num}")
    done

    local rinexver=$(head -1 "$rinexobs" | cut -c 6-6)
    case "$rinexver" in
    "2" )
        local obstypes=$(grep "# / TYPES OF OBSERV" "$rinexobs")
        for sys_num in ${obs_cand[@]}; do
            sys="${sys_num:0:1}"
            sed "1,/END OF HEADER/d" "$rinexobs" | grep -q "$sys[0-9][0-9]" || continue
            obs_cand=("${obs_cand[@]/"${sys}XX"}")
            num="${sys_num:1:2}"
            echo "$obstypes" | grep -q "$num" && obs_cand=("${obs_cand[@]/$sys_num}")
            case "$sys_num" in
            "GC1" ) num="P1" ;;
            "GC2" ) num="P2" ;;
            "RC1" ) num="P1" ;;
            "RC2" ) num="P2" ;;
              *   ) continue ;;
            esac
            echo "$obstypes" | grep -q "$num" && obs_cand=("${obs_cand[@]/$sys_num}")
        done
        ;;
    "3" | "4" )
        local obstypes=$(grep "SYS / # / OBS TYPES" "$rinexobs")
        while IFS= read -r line; do
            if [[ " " != $(echo "$line" | cut -c 1) ]]; then
                sys=$(echo "$line" | cut -c 1) && obs_cand=("${obs_cand[@]/"${sys}XX"}")
            fi
            grep -Eq "^ $sys[0-9][0-9] " "$ctrl_file" || continue
            for sys_num in ${obs_cand[@]}; do
                [ "${sys_num:0:1}" == "$sys" ]        || continue
                obs="${sys_num:1}"
                echo "${freq_cmb[@]}" | grep -Eq "C17" && [ "$sys" == "C" ] && obs="$obs[DPXZ]"
                echo "$line" | grep -Eq "$obs" && obs_cand=("${obs_cand[@]/$sys_num}")
            done
        done <<< "$obstypes"
        ;;
     *  )
        echo -e "$MSGERR unsupported RINEX version (none of 2, 3, 4): $rinexobs"
        return 1
        ;;
    esac

    local avail_sys=("G" "R" "E" "C" "J")
    for sys in ${avail_sys[@]}; do
        grep -Eq "^ $sys[0-9][0-9] " "$ctrl_file" || continue
        if echo "${obs_cand[@]}" | grep -q "${sys}XX"; then
            echo -e "$MSGWAR no observation for GNSS ($sys)"
            continue
        fi
        sys_num=($(echo "${obs_cand[@]}" | grep -o "$sys\w\w" | cut -c 2-3))
        if [ ${#sys_num[@]} -ne 0 ]; then
            echo -e "$MSGWAR no observation for GNSS ($sys): ${sys_num[@]}"
            continue
        fi
    done

    # Prepare RINEX-NAV
    PrepareRinexNav "$mjd_s" "$mjd_e" "$rinex_dir" "$ctrl_file"
    if [ $? -ne 0 ]; then
        echo -e "$MSGERR PrepareRinexNav failed"
        return 1
    fi

    # RINEX-NAV check
    local rinexnav="$rinex_dir/brdm${doy_s}0.${ymd_s:2:2}p"
    [ ! -f "$rinexnav" ] && echo -e "$MSGWAR $rinexnav doesn't exist" && return 1
    local rinexver=$(head -1 "$rinexnav" | cut -c 6-6)
    if [ "$rinexver" -lt "2" -o "$rinexver" -gt 4 ]; then
        echo -e "$MSGERR unsupported RINEX version (none of 2, 3, 4): $rinexnav"
        return 1
    fi

    # Prepare products
    local product_dir=$(get_ctrl "$ctrl_file" "Product directory" | sed "s/^[ ]*//; s/[ ]*$//; s#^~#$HOME#")
    PrepareProducts "$mjd_s" "$mjd_e" "$product_dir" "$ctrl_file" "$rinexobs" "$AR"
    if [ $? -ne 0 ]; then
        echo -e "$MSGERR PrepareProducts failed"
        return 1
    fi

    # Generate binary sp3
    local sp3=$(get_ctrl "$ctrl_file" "Satellite orbit")
    [ $(echo "$sp3" | wc -w) -gt 1 ] && sp3="mersp3_$ymd_s$doy_s"
    cmd="sp3orb $sp3 -cfg $ctrl_file"
    ExecuteWithoutOutput "$cmd"
    if [ $? -ne 0 ]; then
        echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} execution failed"
        return 1
    fi

    # Truncate at midnight (default as NO)
    local tct_opt=$(get_ctrl "$ctrl_file" "Truncate at midnight" | tr 'a-z' 'A-Z')
    if [[ "$tct_opt" == "DEFAULT" ]]; then
        head -1 "$sp3" | grep -q "289   u+U IGS.. FIT  WHU" && tct_opt="NO" || tct_opt="YES"
    fi
    sedi "/^Truncate at midnight/s/ = .*/ = $tct_opt/" "$ctrl_file"

    # Process single site
    ProcessSingleSite "$rinexobs" "$rinexnav" "$ctrl_file" "$mjd_s" "$hms_s" "$mjd_e" "$hms_e" "$site" "$AR"
    if [ $? -ne 0 ]; then
        echo -e "$MSGERR ProcessSingleSession: processing from $ymd_s $doy_s to $ymd_e $doy_e $site failed"
        [ ${DEBUG} == "NO" ] && CleanMid "$ymd_s" "$doy_s"
    else
        echo -e "$MSGSTA ProcessSingleSession from $ymd_s $doy_s to $ymd_e $doy_e done"
    fi
}

ProcessSingleSite() { # purpose : process data of single site
                      # usage   : ProcessSingleSite rinexobs rinexnav config mjd_s hms_s mjd_e hms_e site AR(A/Y/N)
    local rinexobs="$1"
    local rinexnav="$2"
    local config="$3"
    local mjd_s="$4"
    local hms_s="$5"
    local mjd_e="$6"
    local hms_e="$7"
    local site="$8"
    local AR="$9"

    local ydoy_s=($(mjd2ydoy $mjd_s))
    local ydoy_e=($(mjd2ydoy $mjd_e))
    local ymd_s=($(ydoy2ymd ${ydoy_s[*]}))
    local ymd_e=($(ydoy2ymd ${ydoy_s[*]}))

    local year=${ydoy_s[0]}
    local doy=${ydoy_s[1]}
    local ymd=(${ymd_s[@]})

    local freq_cmb=$(get_ctrl "$ctrl_file" "Frequency combination")
    local positioning_mode=$(grep "^ $site [A-Z]" "$config" | awk '{print $2}') # S/P/K/F/L
    local cutoff_elevation=$(grep "^ $site [A-Z]" "$config" | awk '{print $6}') # int, degree

    echo -e "$MSGSTA ProcessSingleSite ${site} from ${ydoy_s[@]} to ${ydoy_e[@]} ..."

    # Compute a priori positions
    echo -e "$MSGSTA Prepare initial position ${site} ..."
    local interval=$(get_ctrl "$config" "Interval")
    ComputeInitialPos "$rinexobs" "$rinexnav" "$mjd_s" "$hms_s" "$mjd_e" "$hms_e" "$site" "$interval" "$positioning_mode"
    if grep -q " NaN " tmp_ComputeInitialPos; then
        echo -e "$MSGERR invalid initial position or sigma: $site"
        return 1
    fi
    local session_time=($(awk '/^Duration/{print $3,$4,$5,$6,$7,$8,$9}' tmp_ComputeInitialPos))
    if [ "$positioning_mode" == "F" ]; then
        if [ -e sit.xyz ] && grep -q "$site" sit.xyz; then
            local initial_pos=($(awk '/'${site}'/{print($2,$3,$4,$5,$6,$7)}' sit.xyz))
            printf " %s%16.4f%16.4f%16.4f%10.6f%10.6f%10.6f\n" ${site} ${initial_pos[*]} > sit.xyz
        else
            local initial_pos=($(snx2sit $site $mjd_s))
            if [ ${#initial_pos[@]} -ne 6 ]; then
                echo -e "$MSGERR ProcessSingleDay: no position or sigma in station coordinate product: $site snx_${year}${doy}"
                return 1
            fi
            printf " %s%16.4f%16.4f%16.4f%10.6f%10.6f%10.6f\n" ${site} ${initial_pos[*]} > sit.xyz
        fi
    else
        awk -v sit=$site '/^Position/{printf(" %s%16.4f%16.4f%16.4f\n",sit,$3,$4,$5)}' tmp_ComputeInitialPos > sit.xyz
        if [ "$positioning_mode" == "L" ]; then
            local product_dir=$(get_ctrl "$ctrl_file" "Product directory" | sed "s/^[ ]*//; s/[ ]*$//; s#^~#$HOME#")
            local kin_pos="kin_${ydoy_s[0]}${ydoy_s[1]}_${site}"
            local pso_pos="pso_${ydoy_s[0]}${ydoy_s[1]}_${site}"
            if [ -e "$product_dir/leo/$pso_pos" ]; then
                mv "${kin_pos}" "${kin_pos}_spp"
                pso2kin.py "$product_dir/leo" "$site" "$mjd_s" "$mjd_e"
                if [ $? -ne 0 -o ! -e "$kin_pos" ]; then
                    echo -e "$MSGWAR ProcessSingleDay: failed to convert PSO to initial coordinates, use results from SPP instead"
                    mv -f "${kin_pos}_spp" "${kin_pos}"
                else
                    rm -f "${kin_pos}_spp"
                fi
            fi
        fi
    fi
    rm -f tmp_ComputeInitialPos
    echo -e "$MSGSTA Prepare initial position ${site} done"

    # Check priori positions
    local xyz=($(awk -v sit=$site '{if($1==sit){print $2,$3,$4}}' sit.xyz))
    if [ -n "$xyz" -a "$positioning_mode" != "L" ]; then
        local blh=($(xyz2blh "${xyz[@]}"))
        if [[ $(echo "${blh[2]} <= -4000" | bc) -eq 1 ]] || \
           [[ $(echo "${blh[2]} >= 20000" | bc) -eq 1 ]]; then
            local blh=$(echo "${blh[2]}/1000" | bc)
            echo -e "$MSGERR ProcessSingleSite: invalid site elevation (out of range from -4 km to +20 km): $site $blh km"
            return 1
        fi
    fi

    # Fill in session time
    if [ ${#session_time[@]} -ne 7 ]; then
        echo -e "$MSGERR ProcessSingleSite: no session time"
        return 1
    else
        session_time="${session_time[@]}"
        sedi "/^Session time/s/ = .*/ = $session_time/" "$ctrl_file"
    fi

    # Editing mode
    local editing_mode=$(get_ctrl "$config" "Strict editing")
    if [[ "$editing_mode" != "YES" ]] && [[ "$editing_mode" != "NO" ]]; then
        echo -e "$MSGERR ProcessSingleSite: illegal editing mode: $editing"
        return 1
    fi
   
    # Tighter threshold for kinematic positioning mode (loose editing mode)
    if [ "$positioning_mode" == "P" -o "$positioning_mode" == "K" ]; then
        local tth_opt="no"
        local lcc_opt="no"
        [[ "$editing_mode" == "NO" ]] && tth_opt="yes"
        [[ "$editing_mode" == "NO" ]] && lcc_opt="lm"
    fi

    # Truncate at midnight
    local tct_opt=$(get_ctrl "$config" "Truncate at midnight" | cut -c 1 | tr 'A-Z' 'a-z')

    # Data preprocess
    echo -e "$MSGSTA Data pre-processing ..."
    local session=$(grep "^Session time" "${config}" | awk '{print $10}')
    local hms=($(grep "^Session time" "${config}" | awk '{print $7,$8,$9}'))
    local rhd_file="log_${year}${doy}_${site}"
    xyz=($(awk -v sit=$site '{if($1==sit){print $2,$3,$4}}' sit.xyz))
    local cmd=""
    if [ "$positioning_mode" == "S" -o "$positioning_mode" == "F" ]; then
        cmd="tedit \"${rinexobs}\" -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
             -xyz ${xyz[*]} -short 1200 -lc_check only -rhd ${rhd_file} -pc_check 300 \
             -elev ${cutoff_elevation} -rnxn \"${rinexnav}\" -freq ${freq_cmb} -trunc_dbd ${tct_opt} -tighter_thre no"
        if [ $mjd_s -le 51666 ]; then
            cmd="tedit \"${rinexobs}\" -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
                 -xyz ${xyz[*]} -short 1200 -lc_check no -rhd ${rhd_file} -pc_check 0 \
                 -elev ${cutoff_elevation} -rnxn \"${rinexnav}\" -freq ${freq_cmb} -trunc_dbd ${tct_opt} -tighter_thre no"
        fi
    elif [ "$positioning_mode" == "P" -o "$positioning_mode" == "K" ]; then
        cmd="tedit \"${rinexobs}\" -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
             -xyz kin_${year}${doy}_${site} -short 120 -lc_check ${lcc_opt} \
             -elev ${cutoff_elevation} -rhd ${rhd_file} -rnxn \"${rinexnav}\" -freq ${freq_cmb} -trunc_dbd ${tct_opt} -tighter_thre ${tth_opt}"
        if [ $mjd_s -le 51666 ]; then
            cmd="tedit \"${rinexobs}\" -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
                 -xyz kin_${year}${doy}_${site} -short 120 -lc_check no \
                 -pc_check 0 -elev ${cutoff_elevation} -rhd ${rhd_file} \
                 -rnxn \"${rinexnav}\" -freq ${freq_cmb} -trunc_dbd ${tct_opt} -tighter_thre ${tth_opt}"
        fi
    elif [ "$positioning_mode" == "L" ]; then
        cmd="tedit \"${rinexobs}\" -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
            -xyz kin_${year}${doy}_${site} -short 120 -lc_check lm \
            -elev ${cutoff_elevation} -rhd ${rhd_file} -rnxn \"${rinexnav}\" -freq ${freq_cmb} -trunc_dbd ${tct_opt} -tighter_thre no"
        if [ $mjd_s -le 51666 ]; then
           cmd="tedit \"${rinexobs}\" -time ${ymd[*]} ${hms[*]} -len ${session} -int ${interval} \
                -xyz kin_${year}${doy}_${site} -short 120 -lc_check lm \
                -pc_check 0 -elev ${cutoff_elevation} -rhd ${rhd_file} \
                -rnxn \"${rinexnav}\" -freq ${freq_cmb} -trunc_dbd ${tct_opt} -tighter_thre no"
        fi
    else
        echo -e "$MSGERR ProcessSingleSite: illegal positioning mode: $positioning_mode"
        return 1
    fi
    cmd=$(tr -s " " <<< "$cmd")
    if [ $(echo "$session > 432000" | bc) -eq 1 ]; then
        Execute "$cmd" || return 1
    else
        ExecuteWithoutOutput "$cmd" || return 1
    fi
    echo -e "$MSGSTA Data pre-processing done"

    # Data clean (iteration)
    echo -e "$MSGSTA Data cleaning ..."
    if [ "$editing_mode" == "YES" ]; then
        local short=$(echo $interval | awk '{printf("%.0f\n", 600/$1)}')
        local jumps=(400 200 100 50)
        local jump_end=50
        if [ "$positioning_mode" == "L" ]; then
            short=$(echo $interval | awk '{printf("%.0f\n", 300/$1)}')
            jumps=(400 200 100 80 60 40)
            jump_end=40
        fi
    else
        local short=0
        local jumps=(400 200 100)
        local jump_end=100
    fi
    local new_rem=100
    local new_amb=100
    for jump in ${jumps[*]}; do
        if [ $new_rem != 0 -o $new_amb != 0 ]; then
          cmd="lsq ${config} \"${rinexobs}\""
          ExecuteWithoutOutput "$cmd" || return 1
        fi
        cmd="redig res_${year}${doy} -jmp $jump -sht $short"
        [ "$positioning_mode" == "L" ] && cmd="$cmd -pce"
        local time=`date +'%Y-%m-%d %H:%M:%S'`
        $cmd > tempout 2>&1
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} execution ok"
        else
            echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} execution failed"
            return 1
        fi
        awk '/%%%\+RMS OF RESIDUALS---PHASE\(MM\)/,/%%%\-RMS OF RESIDUALS---PHASE\(MM\)/{print}' tempout
        new_rem=`awk '/NEWLY REMOVED:/{print $3}' tempout`
        new_amb=`awk '/NEWLY AMBIGUT:/{print $3}' tempout`
        if [ $new_rem == '' -o $new_amb == '' ]; then
            echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} execution failed"
            return 1
        fi
        awk '/NEWLY REMOVED:/{printf "\033[1;34mNewly removed observations\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
        awk '/NEWLY AMBIGUT:/{printf "\033[1;34mNewly inserted ambiguities\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
    done
    niter=0
    while [ $new_rem != 0 -o $new_amb != 0 ]; do
        [ $niter -gt 100 ] && break
        ((niter=niter+1))
        cmd="lsq ${config} \"${rinexobs}\""
        ExecuteWithoutOutput "$cmd" || return 1
        cmd="redig res_${year}${doy} -jmp $jump_end -sht $short"
        [ "$positioning_mode" == "L" ] && cmd="$cmd -pce"
        local time=`date +'%Y-%m-%d %H:%M:%S'`
        $cmd > tempout 2>&1
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} execution ok"
        else
            echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} execution failed"
            return 1
        fi
        new_rem=`awk '/NEWLY REMOVED:/{print $3}' tempout`
        new_amb=`awk '/NEWLY AMBIGUT:/{print $3}' tempout`
        if [ $new_rem == '' -o $new_amb == '' ]; then
            echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} execution failed"
            return 1
        fi
        awk '/%%%\+RMS OF RESIDUALS---PHASE\(MM\)/,/%%%\-RMS OF RESIDUALS---PHASE\(MM\)/{print}' tempout
        awk '/NEWLY REMOVED:/{printf "\033[1;34mNewly removed observations\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
        awk '/NEWLY AMBIGUT:/{printf "\033[1;34mNewly inserted ambiguities\033[0m: %10d%8.2f%%\n",$3,$4}' tempout
    done
    rm -f tempout
    echo -e "$MSGSTA Data cleaning done"

    local vbs_opt=$(get_ctrl "$ctrl_file" "Verbose output" | tr 'a-z' 'A-Z')

    # Ambiguity fixing
    if [ "$AR" != "N" ] && [ -f ?(mer)fcb_${year}${doy} ]; then
        if [ $(grep "# OF AMB RESOLVABLE SAT" amb_${year}${doy} | awk '{print($1)}') -ne 0 ]; then
            cmd="arsig ${config}"
            Execute "$cmd" || return 1
            cmd="lsq ${config} \"${rinexobs}\""
            if [ "$vbs_opt" == "YES" ]; then
                Execute "$cmd" || return 1
            else
                ExecuteWithoutOutput "$cmd" || return 1
            fi
        else
            if [ "$AR" == "Y" ]; then
                echo -e "$MSGERR no resolvable ambiguities"
                echo -e "$MSGINF please check if all required observations and OSBs exist"
                return 1
            elif [ "$AR" == "A" ]; then
                echo -e "$MSGWAR no resolvable ambiguities"
                echo -e "$MSGINF please check if all required observations and OSBs exist"
            fi
        fi
    fi

    echo -e "$MSGSTA Final processing done"

    # Rename result files
    local fn typ types=(rck ztd htg amb res stt cst)
    for typ in ${types[*]}; do
        fn=${typ}_${year}${doy}
        [ -f ${fn} ] && mv -f ${fn} ${fn}_${site}
    done

    echo -e "$MSGSTA ProcessSingleSite ${site} from ${ydoy_s[@]} to ${ydoy_e[@]} done"
}

ComputeInitialPos() { # purpose : compute intial postion with spp
                      # usage   : ComputeInitialPos rinexobs rinexnav mjd_s hms_start mjd_e hms_end site interval mode
    local rinexobs="$1"
    local rinexnav="$2"
    local mjd_s="$3"
    local hms_s="$4"
    local mjd_e="$5"
    local hms_e="$6"
    local site="$7"
    local interval="$8"
    local mode="$9"

    local ydoy_s=($(mjd2ydoy ${mjd_s}))
    local ymd_s=($(ydoy2ymd ${ydoy_s[*]}))
    local ydoy_e=($(mjd2ydoy ${mjd_e}))
    local ymd_e=($(ydoy2ymd ${ydoy_e[*]}))

    local ts="${ymd_s[0]}/${ymd_s[1]}/${ymd_s[2]} $hms_s"
    local te="${ymd_e[0]}/${ymd_e[1]}/${ymd_e[2]} $hms_e"

    local cmd=""
    if [ "$mode" == "K" -o "$mode" == "P" ]; then
        cmd="spp -elev 10 -trop saas -ts $ts -te $te -ti $interval -o kin_${ydoy_s[0]}${ydoy_s[1]}_${site} \"$rinexobs\" \"$rinexnav\""
    elif [ "$mode" == "L" ]; then
        cmd="spp -elev  0 -trop non  -ts $ts -te $te -ti $interval -o kin_${ydoy_s[0]}${ydoy_s[1]}_${site} \"$rinexobs\" \"$rinexnav\""
    else
        cmd="spp -elev 10 -trop saas -ts $ts -te $te -ti $interval \"$rinexobs\" \"$rinexnav\""
    fi

    Execute "$cmd" tmp_ComputeInitialPos || return 1
}

PrepareTables() { # purpose: prepare PRIDE-PPPAR needed tables in working directory
                  # usage  : PrepareTables mjd_s mjd_e table_dir
    local mjd_s="$1"
    local mjd_e="$2"
    local table_dir="$3"

    echo -e "$MSGSTA PrepareTables ..."

    if ! ls $table_dir &>/dev/null; then
        echo -e "$MSGERR PrepareTables: table directory doesn't exist: $table_dir"
        [ "$table_dir" == "Default" ] && echo -e "$MSGINF please specify the table directory in the configuration file"
        return 1
    fi

    local tables=(file_name oceanload orography_ell orography_ell_1x1 gpt3_1.grd)
    for table in ${tables[*]}; do
        if [ ! -f "$table_dir/$table" ]; then
            echo -e "$MSGERR PrepareTables: no such file: $table_dir/$table"
            return 1
        fi
        ln -sf "$table_dir/$table" ./
    done

    # Check leap.sec
    local leapsec="leap.sec"
    local leapsec_ftp="0"
    local leapsec_exi="0"
    local leapsec_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/table/$leapsec"
    if [ -f "$leapsec" ]; then
        sed -n "1p;q" "$leapsec" | grep -q "*" || leapsec_ftp="1"
        grep -q "\-leap sec" "$leapsec"        || leapsec_ftp="1"
    else
        leapsec_exi=$?
    fi
    if [ "$leapsec_ftp" != 0 -o "$leapsec_exi" != 0 ]; then
        rm -f "$leapsec"
        WgetDownload "$leapsec_url"
        if [ ! -f "$leapsec" ]; then
            cp -f "$table_dir/$leapsec" .
        else
            if ! grep -q "\-leap sec" "$leapsec"; then
                echo -e "$MSGWAR PrepareTables: failed to download $leapsec, use $table_dir/$leapsec instead"
                cp -f "$table_dir/$leapsec" .
            fi
            local diff=$(diff "$leapsec" "$table_dir/$leapsec")
            [[ -n "$diff" ]] && cp -f "$leapsec" "$table_dir/"
        fi
    fi
    if ! grep -q "\-leap sec" "$leapsec"; then
        echo -e "$MSGERR PrepareTables: $leapsec invalid or not found"
        echo -e "$MSGINF please download from $leapsec_url to $table_dir for processing"
        return 1
    fi

    # Check sat_parameters
    local satpara="sat_parameters"
    local satpara_ftp="0"
    local satpara_exi="0"
    local satpara_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/table/$satpara"
    if [ -f "$satpara" ]; then
        local tmpymd=($(sed -n "1p;q" "$satpara" | awk '{print(substr($0,56,4),substr($0,61,2),substr($0,64,2))}'))
        local tmpmjd=$(ymd2mjd ${tmpymd[@]})
        [ $? -ne 0 -o "$tmpmjd" -lt "$mjd_e" ] && satpara_ftp="1"
        grep -q "\-prn_indexed" "$satpara"     || satpara_ftp="1"
    else
        satpara_exi=$?
    fi
    if [ "$satpara_ftp" != 0 -o "$satpara_exi" != 0 ]; then
        rm -f "$satpara"
        WgetDownload "$satpara_url"
        if [ ! -f "$satpara" ]; then
            cp -f "$table_dir/$satpara" .
        else
            if ! grep -q "\-prn_indexed" "$satpara"; then
                echo -e "$MSGWAR PrepareTables: failed to download $satpara, use $table_dir/$satpara instead"
                cp -f "$table_dir/$satpara" .
            fi
            local diff=$(diff "$satpara" "$table_dir/$satpara")
            [[ -n "$diff" ]] && cp -f "$satpara" "$table_dir/"
        fi
    fi
    if ! grep -q "\-prn_indexed" "$satpara"; then
        echo -e "$MSGERR PrepareTables: $satpara invalid or not found"
        echo -e "$MSGINF please download from $satpara_url to $table_dir for processing"
        return 1
    else
        local tmpymd=($(sed -n "1p;q" "$satpara" | awk '{print(substr($0,56,4),substr($0,61,2),substr($0,64,2))}'))
        local tmpmjd=$(ymd2mjd ${tmpymd[@]})
        if [ $? -ne 0 -o "$tmpmjd" -lt "$mjd_e" ]; then
            echo -e "$MSGWAR PrepareTables: outdated $satpara"
            echo -e "$MSGINF please update from $satpara_url to $table_dir for processing"
        fi
    fi

    echo -e "$MSGSTA PrepareTables done"
}

PrepareRinexNav() { # purpose : prepare RINEX multi-systems broadcast ephemerides
                    # usage   : PrepareRinexNav mjd_s mjd_e rinex_dir config
    local mjd_s="$1"
    local mjd_e="$2"
    local rinex_dir="$3"
    local config="$4"

    echo -e "$MSGSTA PrepareRinexNav ..."

    for mjd in $(seq $mjd_s $mjd_e); do
        local ydoy=($(mjd2ydoy $mjd))
        local year=${ydoy[0]}
        local doy=${ydoy[1]}
        local rinexnav="brdm${doy}0.${year:2:2}p"

        # Try downloading hourly navigation file when processing current day's data
        if [ $(date -u +"%Y%j") -eq "$year$doy" ]; then
            local navgps="hour${doy}0.${year:2:2}n" && rm -f "$navgps"
            local urlnav="ftp://igs.gnsswhu.cn/pub/gps/data/hourly/${year}/${doy}/${navgps}.gz"
            WgetDownload "$urlnav"
            if [ $? -eq 0 ]; then
                gunzip -f ${navgps}.gz
            else
                echo -e "$MSGWAR failed to download hourly GPS navigation file: $navgps"
            fi
            local navglo="hour${doy}0.${year:2:2}g" && rm -f "$navglo"
            local urlnav="ftp://igs.gnsswhu.cn/pub/gps/data/hourly/${year}/${doy}/${navglo}.gz"
            WgetDownload "$urlnav"
            if [ $? -eq 0 ]; then
                gunzip -f ${navglo}.gz
            else
                echo -e "$MSGWAR failed to download hourly GLONASS navigation file: $navglo"
            fi
            if [ -f "$navgps" -a -f "$navglo" ]; then
                echo -e "$MSGSTA Merging $rinexnav ..."
                merge2brdm.py "$navgps" "$navglo" && mv -f "$rinexnav" "$rinex_dir"
                if [ $? -ne 0 -o ! -f "$rinex_dir/$rinexnav" ]; then
                    echo -e "$MSGERR failed to merge hourly RINEX navigation files: $navgps $navglo -> $rinexnav"
                    return 1
                fi
            else
                [ -f "$navglo" ] && mv -f "$navglo" "$rinex_dir/$rinexnav"
                [ -f "$navgps" ] && mv -f "$navgps" "$rinex_dir/$rinexnav"
                if [ ! -f "$rinex_dir/$rinexnav" ]; then
                    echo -e "$MSGERR failed to download hourly RINEX navigation file: $rinex_dir/$rinexnav"
                    echo -e "$MSGINF please download from $urlnav to $rinex_dir for processing"
                    return 1
                fi
            fi
        fi

        # Try finding brdm from RINEX directory
        if [ ! -f "$rinex_dir/$rinexnav" ]; then
            local tmpnav="BRDC00IGS_R_${year}${doy}0000_01D_MN.rnx"
            [ -f "$rinex_dir/$tmpnav" -a ! -f "$rinex_dir/$rinexnav" ] && mv -f "$rinex_dir/$tmpnav" "$rinex_dir/$rinexnav"
            local tmpnav="BRDC00IGN_R_${year}${doy}0000_01D_MN.rnx"
            [ -f "$rinex_dir/$tmpnav" -a ! -f "$rinex_dir/$rinexnav" ] && mv -f "$rinex_dir/$tmpnav" "$rinex_dir/$rinexnav"
            local tmpnav="BRDM00DLR_S_${year}${doy}0000_01D_MN.rnx"
            [ -f "$rinex_dir/$tmpnav" -a ! -f "$rinex_dir/$rinexnav" ] && mv -f "$rinex_dir/$tmpnav" "$rinex_dir/$rinexnav"
        fi

        # Try downloading brdm
        if [ ! -f "$rinex_dir/$rinexnav" ]; then
            if [ $year -ge 2016 ]; then
                local tmpnav="BRDC00IGS_R_${year}${doy}0000_01D_MN.rnx"
                local urlnav="ftp://igs.gnsswhu.cn/pub/gps/data/daily/${year}/${doy}/${year:2:2}p/${tmpnav}.gz"
                WgetDownload "$urlnav"
                if [ $? -eq 0 ]; then
                    gunzip -f ${tmpnav}.gz && mv -f "$tmpnav" "$rinex_dir/$rinexnav"
                else
                    tmpnav="BRDC00IGN_R_${year}${doy}0000_01D_MN.rnx"
                    urlnav="ftp://igs.ign.fr/pub/igs/data/${year}/${doy}/${tmpnav}.gz"
                    WgetDownload "$urlnav"
                    if [ $? -eq 0 ]; then
                        gunzip -f ${tmpnav}.gz && mv -f "$tmpnav" "$rinex_dir/$rinexnav"
                    else
                        echo -e "$MSGWAR failed to download RINEX navigation file: $tmpnav"
                    fi
                fi
            fi
        fi

        # Try downloading GPS and GLONASS brdc
        if [ ! -f "$rinex_dir/$rinexnav" ]; then
            local navgps="brdc${doy}0.${year:2:2}n"
            if [ ! -f "$navgps" ]; then
                local urlnav="ftp://igs.gnsswhu.cn/pub/gps/data/daily/${year}/${doy}/${year:2:2}n/${navgps}.Z"
                CopyOrDownloadProduct "$rinex_dir/$navgps"
                if [ $? -ne 0 ]; then
                    WgetDownload "$urlnav" && gunzip -f "${navgps}.Z"
                    if [ $? -ne 0 ]; then
                        echo -e "$MSGWAR failed to download GPS navigation file: $navgps"
                    fi
                fi
            fi
            local navglo="brdc${doy}0.${year:2:2}g"
            if [ ! -f "$navglo" ]; then
                local urlnav="ftp://igs.gnsswhu.cn/pub/gps/data/daily/${year}/${doy}/${year:2:2}g/${navglo}.Z"
                CopyOrDownloadProduct "$rinex_dir/$navglo"
                if [ $? -ne 0 ]; then
                    WgetDownload "$urlnav" && gunzip -f "${navglo}.Z"
                    if [ $? -ne 0 ]; then
                        echo -e "$MSGWAR failed to download GLONASS navigation file: $navglo"
                    fi
                fi
            fi
            if [ -f "$navgps" -a -f "$navglo" ]; then
                merge2brdm.py "$navgps" "$navglo" && mv -f "$rinexnav" "$rinex_dir"
                if [ $? -ne 0 -o ! -f "$rinex_dir/$rinexnav" ]; then
                    echo -e "$MSGERR failed to merge RINEX navigation files: $navgps $navglo -> $rinexnav"
                    return 1
                fi
            else
                [ -f "$navglo" ] && mv -f "$navglo" "$rinex_dir/$rinexnav"
                [ -f "$navgps" ] && mv -f "$navgps" "$rinex_dir/$rinexnav"
                if [ ! -f "$rinex_dir/$rinexnav" ]; then
                    echo -e "$MSGERR failed to download RINEX navigation file: $rinex_dir/$rinexnav"
                    echo -e "$MSGINF please download from $urlnav to $rinex_dir for processing"
                    return 1
                fi
            fi
        fi

        # Check brdm for each GNSS
        local sys="G" nsys="0"
        local rinexver=$(head -1 "$rinex_dir/$rinexnav" | cut -c 6-6)
        case "$rinexver" in
        "2" )
            head -1 "$rinex_dir/$rinexnav" | grep -Eq "GLO"       && sys="R"
            head -1 "$rinex_dir/$rinexnav" | grep -Eq "GAL"       && sys="E"
            head -1 "$rinex_dir/$rinexnav" | grep -Eq "(COM|BEI)" && sys="C"
            head -1 "$rinex_dir/$rinexnav" | grep -Eq "QZS"       && sys="J"
            grep -q "^ $sys[0-9][0-9] " "$config" && nsys=$[$nsys+1]
            echo -e "$MSGWAR using single-GNSS ($sys) RINEX navigation file: $rinexnav"
            ;;
        "3" | "4" )
            local avail_sys=("G" "R" "E" "C" "J")
            for sys in ${avail_sys[@]}; do
               grep -Eq "^ $sys[0-9][0-9] " "$config" || continue
               grep -Eq "^$sys[ 0-9][0-9] " "$rinex_dir/$rinexnav" && nsys=$[$nsys+1]
               if [ $? -ne 0 ]; then
                   echo -e "$MSGWAR no navigation message for GNSS ($sys): $rinexnav"
               fi
            done
            ;;
        esac
        if [ "$nsys" -eq 0 ]; then
            echo -e "$MSGERR all satellites in RINEX navigation file are disabled: $rinexnav"
            exit 1
        fi
    done

    echo -e "$MSGSTA PrepareRinexNav done"
}

PrepareProducts() { # purpose : prepare PRIDE-PPPAR needed products in working directory
                    # usage   : PrepareProducts mjd_s mjd_e product_dir config rinexobs AR(A/Y/N)
    local mjd_s="$1"
    local mjd_e="$2"
    local product_dir="$3"
    local config="$4"
    local rinexobs="$5"
    local AR="$6"

    local ydoy_s=($(mjd2ydoy $mjd_s))
    local ymd_s=($(ydoy2ymd ${ydoy_s[*]}))
    local doy_s=${ydoy_s[1]}

    local ign_priority_path="$(dirname $(which pdp3))/.ign_priority"

    echo -e "$MSGSTA PrepareProducts ..."

    if [ "$product_dir" == "Default" ]; then
        product_dir="$(rdlk ..)/product/"
        sedi "/^Product directory/s/ = .*/ = ${product_dir//\//\\/}/" "$config"
    fi

    local product_cmn_dir="$product_dir/common"
    local product_ion_dir="$product_dir/ion"
    local product_vmf_dir="$product_dir/vmf"
    local product_ssc_dir="$product_dir/ssc"
    local product_leo_dir="$product_dir/leo"

    [ "$OFFLINE" == "NO" ] && mkdir -p "$product_cmn_dir"

    # Satellite orbit
    local custom_pro_sp3=$(get_ctrl "$config" "Satellite orbit")
    if [ "$custom_pro_sp3" != Default ]; then
        local sp3="$custom_pro_sp3"
        for sp3 in $custom_pro_sp3; do
            CopyOrDownloadProduct "$product_cmn_dir/$sp3"
            if [ $? -ne 0 ]; then
                CopyOrDownloadProduct "$product_dir/$sp3"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: no satellite orbit product: $sp3"
                    return 1
                fi
            fi          
        done
        local argnum=$(echo "$custom_pro_sp3" | wc -w)
        if [ $argnum -gt 1 ]; then    
            sp3="mersp3_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_sp3" "$sp3"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge satellite orbit products: $custom_pro_sp3 -> $sp3"
                return 1
            fi
            for tmp in $(echo "$custom_pro_sp3"); do
                rm -f "$tmp"
            done
        fi
    else
        local custom_pro_sp3=""
        for mjd in $(seq $mjd_s $mjd_e); do
            if [ "$mjd_s" -lt 49718 ]; then
                echo -e "$MSGERR no available satellite product for dates before MJD 49718"
                return 1
            fi
            local ydoy=($(mjd2ydoy $mjd))
            local wkdow=($(mjd2wkdow $mjd_s))
            local urls=(
                "ftp://igs.ign.fr/pub/igs/products/mgex/${wkdow[0]}/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_05M_ORB.SP3.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_05M_ORB.SP3.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_01M_ORB.SP3.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/IGS2R03FIN_${ydoy[0]}${ydoy[1]}0000_01D_05M_ORB.SP3.gz"
            )
            for url in ${urls[@]}; do
                if [[ "$url" =~ igs.ign.fr ]]; then
                    [ -e "$ign_priority_path" ] && [ "${wkdow[0]}" -ge 2290 ] || continue
                fi
                local cmp=$(basename "$url")
                local sp3="${cmp/\.[gZ]*/}"
                CopyOrDownloadProduct "$product_cmn_dir/$sp3"        && break
                CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url" && break
            done
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$sp3" ]; then
                local mjd_t=$(ymd2mjd $(date +"%Y %m %d"))
                if [ $mjd_e -gt $(($mjd_t-3)) ] && [ "$USERTS" == "YES" ]; then
                    local url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/WUM0MGXRTS_${ydoy[0]}${ydoy[1]}0000_01D_05M_ORB.SP3.gz"
                    local cmp=$(basename "$url")
                    local sp3="${cmp/\.[gZ]*/}"
                    echo -e "$MSGWAR PrepareProducts: failed to download RAP satellite orbit product $cmp, try downloading RTS products"
                    if [ -f "$product_cmn_dir/$cmp" ]; then
                        size_last=$(ls -l "$product_cmn_dir/$cmp" | awk '{print($5)}')
                        size_next=$(curl "$(dirname $url)/" | grep "$sp3" | awk '{print($5)}')
                        if [ $? -eq 0 ]; then
                            if [ "$size_next" -gt "$size_last" ]; then
                                rm -f "$sp3"* "$product_cmn_dir/$sp3"*
                            fi
                        fi
                    fi
                    CopyOrDownloadProduct "$product_cmn_dir/$sp3"
                    [ $? -ne 0 ] && CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url"
                    ## No satellite attitude product for RTS
                    custom_pro_att="NONE"
                    local att="$custom_pro_att"
                    sedi "/Quaternions/s/Default/$att/g" "$config"
                else
                    [ ${ydoy[0]} -ge 2020 ] && local url="${urls[1]}" || local url="${urls[${#urls[@]}-1]}"
                    local cmp=$(basename "$url")
                    local sp3="${cmp/\.[gZ]*/}"
                fi
            fi
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$sp3" ]; then
                echo -e "$MSGERR PrepareProducts: failed to download satellite orbit product: $cmp"
                echo -e "$MSGINF please download from $url to $product_cmn_dir for processing"
                return 1
            fi
            sedi "/Satellite orbit/s/Default/$sp3 &/" "$config"
            custom_pro_sp3="$custom_pro_sp3 $sp3"
        done
        sedi "/Satellite orbit/s/Default//g" "$config"
        local argnum=$(echo "$custom_pro_sp3" | wc -w)
        if [ $argnum -gt 1 ]; then
            sp3="mersp3_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_sp3" "$sp3"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge satellite orbit products: $custom_pro_sp3 -> $sp3"
                return 1
            fi
            for tmp in $(echo "$custom_pro_sp3"); do
                rm -f "$tmp"
            done
        fi        
    fi

    # Satellite clock
    local custom_pro_clk=$(get_ctrl "$config" "Satellite clock")
    if [ "$custom_pro_clk" != Default ]; then        
        local clk="$custom_pro_clk"
        for clk in $custom_pro_clk; do
            CopyOrDownloadProduct "$product_cmn_dir/$clk"
            if [ $? -ne 0 ]; then
                CopyOrDownloadProduct "$product_dir/$clk"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: no satellite clock product: $clk"
                    return 1
                fi
            fi          
        done
        local argnum=$(echo "$custom_pro_clk" | wc -w)
        if [ $argnum -gt 1 ]; then    
            clk="mersck_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_clk" "$clk"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge satellite clock products: $custom_pro_clk -> $clk"
                return 1
            fi
            for tmp in $(echo "$custom_pro_clk"); do
                rm -f "$tmp"
            done
        fi        
    else
        local custom_pro_clk=""
        for mjd in $(seq $mjd_s $mjd_e); do
            if [ "$mjd_s" -lt 49718 ]; then
                echo -e "$MSGERR no available satellite product for dates before MJD 49718"
                return 1
            fi
            local ydoy=($(mjd2ydoy $mjd))
            local wkdow=($(mjd2wkdow $mjd_s))
            local urls=(
                "ftp://igs.ign.fr/pub/igs/products/mgex/${wkdow[0]}/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_30S_CLK.CLK.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/clock/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_30S_CLK.CLK.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/clock/IGS2R03FIN_${ydoy[0]}${ydoy[1]}0000_01D_30S_CLK.CLK.gz"
            )
            for url in ${urls[@]}; do
                if [[ "$url" =~ igs.ign.fr ]]; then
                    [ -e "$ign_priority_path" ] && [ "${wkdow[0]}" -ge 2290 ] || continue
                fi
                local cmp=$(basename "$url")
                local clk="${cmp/\.[gZ]*/}"
                CopyOrDownloadProduct "$product_cmn_dir/$clk"        && break
                CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url" && break
            done
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$clk" ]; then
                local mjd_t=$(ymd2mjd $(date +"%Y %m %d"))
                if [ $mjd_e -gt $(($mjd_t-3)) ] && [ "$USERTS" == "YES" ]; then
                    local url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/clock/WUM0MGXRTS_${ydoy[0]}${ydoy[1]}0000_01D_05S_CLK.CLK.gz"
                    local cmp=$(basename "$url")
                    local clk="${cmp/\.[gZ]*/}"
                    echo -e "$MSGWAR PrepareProducts: failed to download RAP satellite clock product $cmp, try downloading RTS products"
                    if [ -f "$product_cmn_dir/$cmp" ]; then
                        size_last=$(ls -l "$product_cmn_dir/$cmp" | awk '{print($5)}')
                        size_next=$(curl "$(dirname $url)/" | grep "$clk" | awk '{print($5)}')
                        if [ $? -eq 0 ]; then
                            if [ "$size_next" -gt "$size_last" ]; then
                                rm -f "$clk"* "$product_cmn_dir/$clk"*
                            fi
                        fi
                    fi
                    CopyOrDownloadProduct "$product_cmn_dir/$clk"
                    [ $? -ne 0 ] && CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url"
                else
                    [ ${ydoy[0]} -ge 2020 ] && local url="${urls[1]}" || local url="${urls[${#urls[@]}-1]}"
                    local cmp=$(basename "$url")
                    local clk="${cmp/\.[gZ]*/}"
                fi
            fi
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$clk" ]; then
                echo -e "$MSGERR PrepareProducts: failed to download satellite clock product: $cmp"
                echo -e "$MSGINF please download from $url to $product_cmn_dir for processing"
                return 1
            fi
            sedi "/Satellite clock/s/Default/$clk &/" "$config"
            custom_pro_clk="$custom_pro_clk $clk"
        done
        sedi "/Satellite clock/s/Default//g" "$config"
        local argnum=$(echo "$custom_pro_clk" | wc -w)
        if [ $argnum -gt 1 ]; then    
            clk="mersck_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_clk" "$clk"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge satellite clock products: $custom_pro_clk -> $clk"
                return 1
            fi
            for tmp in $(echo "$custom_pro_clk"); do
                rm -f "$tmp"
            done
        fi
    fi

    # Earth rotation parameters
    local custom_pro_erp=$(get_ctrl "$config" "ERP")
    if [ "$custom_pro_erp" != Default ]; then
        local erp="$custom_pro_erp"
        for erp in $custom_pro_erp; do
            CopyOrDownloadProduct "$product_cmn_dir/$erp"
            if [ $? -ne 0 ]; then
                CopyOrDownloadProduct "$product_dir/$erp"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: no ERP product: $erp"
                    return 1
                fi
            fi
        done
        local argnum=$(echo "$custom_pro_erp" | wc -w)
        if [ $argnum -gt 1 ]; then
            erp="mererp_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_erp" "$erp"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge ERP products: $custom_pro_erp -> $erp"
                return 1
            fi
            for tmp in $(echo "$custom_pro_erp"); do
                rm -f "$tmp"
            done
        fi
    else
        local custom_pro_erp=""
        local mjd_seq=($(seq $mjd_s $mjd_e))
        [ "$mode" == "L" ] && mjd_seq=($(seq $[$mjd_s-1] $[$mjd_e+1]))
        for mjd in "${mjd_seq[@]}"; do
            if [ "$mjd_s" -lt 49718 ]; then
                echo -e "$MSGERR no available satellite product for dates before MJD 49718"
                return 1
            fi
            local ydoy=($(mjd2ydoy $mjd))
            local wkdow=($(mjd2wkdow $mjd_s))
            local urls=(
                "ftp://igs.ign.fr/pub/igs/products/mgex/${wkdow[0]}/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_01D_ERP.ERP.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_01D_ERP.ERP.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/COD0R03FIN_${ydoy[0]}${ydoy[1]}0000_01D_01D_ERP.ERP.gz"
            )
            for url in ${urls[@]}; do
                if [[ "$url" =~ igs.ign.fr ]]; then
                    [ -e "$ign_priority_path" ] && [ "${wkdow[0]}" -ge 2290 ] || continue
                fi
                local cmp=$(basename "$url")
                local erp="${cmp/\.[gZ]*/}"
                CopyOrDownloadProduct "$product_cmn_dir/$erp"        && break
                CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url" && break
            done
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$erp" ]; then
                local mjd_t=$(ymd2mjd $(date +"%Y %m %d"))
                if [ $mjd_e -gt $(($mjd_t-3)) ] && [ "$USERTS" == "YES" ]; then
                    local url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/WUM0MGXRTS_${ydoy[0]}${ydoy[1]}0000_01D_01D_ERP.ERP.gz"
                    local cmp=$(basename "$url")
                    local erp="${cmp/\.[gZ]*/}"
                    echo -e "$MSGWAR PrepareProducts: failed to download RAP ERP product $cmp, try downloading RTS products"
                    if [ -f "$product_cmn_dir/$cmp" ]; then
                        size_last=$(ls -l "$product_cmn_dir/$cmp" | awk '{print($5)}')
                        size_next=$(curl "$(dirname $url)/" | grep "$erp" | awk '{print($5)}')
                        if [ $? -eq 0 ]; then
                            if [ "$size_next" -gt "$size_last" ]; then
                                rm -f "$erp"* "$product_cmn_dir/$erp"*
                            fi
                        fi
                    fi
                    CopyOrDownloadProduct "$product_cmn_dir/$erp"
                    [ $? -ne 0 ] && CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url"
                else
                    [ ${ydoy[0]} -ge 2020 ] && local url="${urls[1]}" || local url="${urls[${#urls[@]}-1]}"
                    local cmp=$(basename "$url")
                    local erp="${cmp/\.[gZ]*/}"
                fi
            fi
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$erp" ]; then
                echo -e "$MSGERR PrepareProducts: failed to download ERP product: $cmp"
                echo -e "$MSGINF please download from $url to $product_cmn_dir for processing"
                return 1
            fi
            sedi "/ERP/s/Default/$erp &/" "$config"
            custom_pro_erp="$custom_pro_erp $erp"
        done
        sedi "/ERP/s/Default//g" "$config"
        local argnum=$(echo "$custom_pro_erp" | wc -w)
        if [ $argnum -gt 1 ]; then    
            erp="mererp_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_erp" "$erp"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge ERP products: $custom_pro_erp -> $erp"
                return 1
            fi
            for tmp in $(echo "$custom_pro_erp"); do
                rm -f "$tmp"
            done 
        fi
    fi

    # Quaternions
    local custom_pro_att=$(get_ctrl "$config" "Quaternions")
    if [ "$custom_pro_att" != Default ]; then
        local att="${custom_pro_att}"
        if [ "$(echo "$custom_pro_att" | tr 'a-z' 'A-Z')" != NONE ]; then
            for att in $custom_pro_att; do
                CopyOrDownloadProduct "$product_cmn_dir/$att"
                if [ $? -ne 0 ]; then
                    CopyOrDownloadProduct "$product_dir/$att"
                    if [ $? -ne 0 ]; then
                        echo -e "$MSGERR PrepareProducts: no satellit attitude product: $att"
                        return 1
                    fi
                fi
            done
            local argnum=$(echo "$custom_pro_att" | wc -w)
            if [ $argnum -gt 1 ]; then    
                att="meratt_${ymd_s}${doy_s}"
                MergeFiles "$(pwd)" "$custom_pro_att" "$att"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: failed to merge satellite attitude products: $custom_pro_att -> $att"
                    return 1
                fi
                for tmp in $(echo "$custom_pro_att"); do
                    rm -f "$tmp"
                done
            fi
        fi
    else
        local custom_pro_att=""
        for mjd in $(seq $mjd_s $mjd_e); do
            if [ "$mjd_s" -lt 49718 ]; then
                echo -e "$MSGERR no available satellite product for dates before MJD 49718"
                return 1
            fi
            local ydoy=($(mjd2ydoy $mjd))
            local wkdow=($(mjd2wkdow $mjd_s))
            local urls=(
                "ftp://igs.ign.fr/pub/igs/products/mgex/${wkdow[0]}/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_30S_ATT.OBX.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_30S_ATT.OBX.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/orbit/IGS2R03FIN_${ydoy[0]}${ydoy[1]}0000_01D_30S_ATT.OBX.gz"
            )
            for url in ${urls[@]}; do
                if [[ "$url" =~ igs.ign.fr ]]; then
                    [ -e "$ign_priority_path" ] && [ "${wkdow[0]}" -ge 2290 ] || continue
                fi
                local cmp=$(basename "$url")
                local att="${cmp/\.[gZ]*/}"
                CopyOrDownloadProduct "$product_cmn_dir/$att"        && break
                CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url" && break
            done
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$att" ]; then
                [ ${ydoy[0]} -ge 2020 ] && local url="${urls[1]}" || local url="${urls[${#urls[@]}-1]}"
                local cmp=$(basename "$url")
                local att="${cmp/\.[gZ]*/}"
                echo -e "$MSGWAR PrepareProducts: failed to download satellite attitude product: $cmp"
                echo -e "$MSGINF please download from $url to $product_cmn_dir for processing"
                break
            fi
            sedi "/Quaternions/s/Default/$att &/" "$config"
            custom_pro_att="$custom_pro_att $att"
        done
        sedi "/Quaternions/s/Default//g" "$config"
        local argnum=$(echo "$custom_pro_att" | wc -w)
        if [ $argnum -gt 1 ]; then
            att="meratt_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_att" "$att"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge satellite attitude products: $custom_pro_att -> $att"
                return 1
            fi
            for tmp in $(echo "$custom_pro_att"); do
                rm -f "$tmp"
            done
        elif [ $argnum -eq 0 ]; then
            att="NONE"
            sedi "/Quaternions/s/Default/$att/g" "$config"
        fi
    fi

    # Code/phase bias
    local custom_pro_fcb=$(get_ctrl "$config" "Code/phase bias")
    if [ "$custom_pro_fcb" != Default ]; then
        local fcb="$custom_pro_fcb"
        if [ "$(echo "$custom_pro_fcb" | tr 'a-z' 'A-Z')" != NONE ]; then
            for fcb in $custom_pro_fcb; do
                CopyOrDownloadProduct "$product_cmn_dir/$fcb"
                if [ $? -ne 0 ]; then
                    CopyOrDownloadProduct "$product_dir/$fcb"
                    if [ $? -ne 0 ]; then
                        if [ "$AR" == "Y" ]; then
                            echo -e "$MSGERR PrepareProducts: no satellite code/phase bias product: $fcb"
                            return 1
                        else
                            echo -e "$MSGWAR PrepareProducts: no satellite code/phase bias product: $fcb"
                        fi
                    fi
                fi
            done
            local argnum=$(echo "$custom_pro_fcb" | wc -w)
            if [ $argnum -gt 1 ]; then    
                fcb="merfcb_${ymd_s}${doy_s}"
                MergeFiles "$(pwd)" "$custom_pro_fcb" "$fcb"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: failed to merge satellite code/phase bias products: $custom_pro_fcb -> $fcb"
                    return 1
                fi
                for tmp in $(echo "$custom_pro_fcb"); do
                    rm -f "$tmp"
                done
            fi
        fi
    else
        local custom_pro_fcb=""
        for mjd in $(seq $mjd_s $mjd_e); do
            if [ "$mjd_s" -lt 51544 ]; then
                if [ "$AR" == "Y" ]; then
                    echo -e "$MSGERR no available satellite code/phase bias product for dates before MJD 51544"
                    return 1
                else
                    echo -e "$MSGINF no available satellite code/phase bias product for dates before MJD 51544"
                fi
                custom_pro_fcb="NONE"
                local fcb="$custom_pro_fcb"
                sedi "/Code\/phase bias/s/Default/$fcb/g" "$config"
                break
            fi
            local ydoy=($(mjd2ydoy $mjd))
            local wkdow=($(mjd2wkdow $mjd_s))
            local urls=(
                "ftp://igs.ign.fr/pub/igs/products/mgex/${wkdow[0]}/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_01D_OSB.BIA.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/bias/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_01D_OSB.BIA.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/bias/WUM0MGXRAP_${ydoy[0]}${ydoy[1]}0000_01D_01D_ABS.BIA.gz"
                "ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/bias/IGS2R03FIN_${ydoy[0]}${ydoy[1]}0000_01D_01D_OSB.BIA.gz"
            )
            for url in ${urls[@]}; do
                if [[ "$url" =~ igs.ign.fr ]]; then
                    [ -e "$ign_priority_path" ] && [ "${wkdow[0]}" -ge 2290 ] || continue
                fi
                local cmp=$(basename "$url")
                local fcb="${cmp/\.[gZ]*/}"
                CopyOrDownloadProduct "$product_cmn_dir/$fcb"        && break
                CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url" && break
            done
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$fcb" ]; then
                local mjd_t=$(ymd2mjd $(date +"%Y %m %d"))
                if [ $mjd_e -gt $(($mjd_t-3)) ] && [ "$USERTS" == "YES" ]; then
                    local url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${ydoy[0]}/bias/WUM0MGXRTS_${ydoy[0]}${ydoy[1]}0000_01D_05M_OSB.BIA.gz"
                    local cmp=$(basename "$url")
                    local fcb="${cmp/\.[gZ]*/}"
                    echo -e "$MSGWAR PrepareProducts: failed to download RAP satellite code/phase bias product $cmp, try downloading RTS products"
                    if [ -f "$product_cmn_dir/$cmp" ]; then
                        size_last=$(ls -l "$product_cmn_dir/$cmp" | awk '{print($5)}')
                        size_next=$(curl "$(dirname $url)/" | grep "$fcb" | awk '{print($5)}')
                        if [ $? -eq 0 ]; then
                            if [ "$size_next" -gt "$size_last" ]; then
                                rm -f "$fcb"* "$product_cmn_dir/$fcb"*
                            fi
                        fi
                    fi
                    CopyOrDownloadProduct "$product_cmn_dir/$fcb"
                    [ $? -ne 0 ] && CopyOrDownloadProduct "$product_cmn_dir/$cmp" "$url"
                else
                    [ ${ydoy[0]} -ge 2020 ] && local url="${urls[1]}" || local url="${urls[${#urls[@]}-1]}"
                    local cmp=$(basename "$url")
                    local fcb="${cmp/\.[gZ]*/}"
                fi
            fi
            [ -f "$cmp" ] && gunzip -f "$cmp"
            if [ ! -f "$fcb" ]; then
                if [ "$AR" == "Y" ]; then
                    echo -e "$MSGERR PrepareProducts: failed to download satellite code/phase bias lock product: $cmp"
                    echo -e "$MSGINF please download from $url to $product_cmn_dir for processing"
                    return 1
                else
                    echo -e "$MSGWAR PrepareProducts: failed to download satellite code/phase bias lock product: $cmp"
                    echo -e "$MSGINF please download from $url to $product_cmn_dir for processing"
                fi
                break
            fi
            sedi "/Code\/phase bias/s/Default/$fcb &/" "$config"
            custom_pro_fcb="$custom_pro_fcb $fcb"
        done
        sedi "/Code\/phase bias/s/Default//g" "$config"
        local argnum="$(echo $custom_pro_fcb | wc -w)"
        if [ $argnum -gt 1 ]; then
            fcb="merfcb_${ymd_s}${doy_s}"
            MergeFiles "$(pwd)" "$custom_pro_fcb" "$fcb"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge satellite code/phase bias products: $custom_pro_fcb -> $fcb"
                return 1
            fi
            for tmp in $(echo "$custom_pro_fcb"); do
                rm -f "$tmp"
            done
        elif [ $argnum -eq 0 ]; then
            fcb="NONE"
            sedi "/Code\/phase bias/s/Default/$fcb/g" "$config"
        fi
    fi

    # LEO quaternions
    local custom_pro_lat=$(get_ctrl "$config" "LEO quaternions")
    if [ "$mode" == "L" -a "$(echo "$custom_pro_lat" | tr 'a-z' 'A-Z')" != NONE ]; then
        mkdir -p "$product_leo_dir"
        if [ "$custom_pro_lat" != Default ]; then
            local lat="$custom_pro_lat"
            for lat in $custom_pro_lat; do
                CopyOrDownloadProduct "$product_leo_dir/$lat"
                if [ $? -ne 0 ]; then
                    CopyOrDownloadProduct "$product_dir/$lat"
                    if [ $? -ne 0 ]; then
                        echo -e "$MSGERR PrepareProducts: no LEO attitude product: $lat"
                        return 1
                    fi
                fi
            done
            local argnum=$(echo "$custom_pro_lat" | wc -w)
            if [ $argnum -gt 1 ]; then
                lat="merlat_${ymd_s}${doy_s}"
                MergeFiles "$(pwd)" "$custom_pro_lat" "$lat"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: failed to merge LEO attitude products: $custom_pro_lat -> $lat"
                    return 1
                fi
                for tmp in $(echo "$custom_pro_lat"); do
                    rm -f "$tmp"
                done
            fi
        else
            local custom_pro_lat=""
            for mjd in $(seq $mjd_s $mjd_e); do
                local ydoy=($(mjd2ydoy $mjd))
                lat="lat_${ydoy[0]}${ydoy[1]}_${site}"
                CopyOrDownloadProduct "$product_leo_dir/$lat"
                if [ $? -ne 0 ]; then
                    CopyOrDownloadProduct "$product_dir/$lat"
                    if [ $? -ne 0 ]; then
                        echo -e "$MSGWAR PrepareProducts: no LEO attitude product: $lat"
                        echo -e "$MSGINF please download with prepare_leodata.sh to $product_leo_dir for processing"
                        custom_pro_lat="NONE"
                        local lat="$custom_pro_lat"
                        sedi "/LEO quaternions/s/Default/$lat/" "$config"
                        break
                    fi
                fi
                sedi "/LEO quaternions/s/Default/$lat &/" "$config"
                custom_pro_lat="$custom_pro_lat $lat"
            done
            sedi "/LEO quaternions/s/Default//g" "$config"
            local argnum=$(echo "$custom_pro_lat" | wc -w)
            if [ $argnum -gt 1 ]; then
                lat="merlat_${ymd_s}${doy_s}"
                MergeFiles "$(pwd)" "$custom_pro_lat" "$lat"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: failed to merge LEO attitude products: $custom_pro_lat -> $lat"
                    return 1
                fi
                for tmp in $(echo "$custom_pro_lat"); do
                    rm -f "$tmp"
                done
            elif [ $argnum -eq 0 ]; then
                lat="NONE"
                sedi "/LEO quaternions/s/Default/$lat/g" "$config"
            fi
        fi
    else
        custom_pro_lat="NONE"
        local lat="$custom_pro_lat"
        sedi "/LEO quaternions/s/Default/$lat/" "$config"
    fi

    # Check satellite orbit
    if [ -f "$sp3" ]; then
        head -1 "$sp3" | grep -q "^#a"
        if [ $? -eq 0 ]; then
            echo -e "$MSGERR unsupproted SP3 version (#a): $custom_pro_sp3"
            return 1
        fi
        local avail_sys=("G" "R" "E" "C" "J")
        for sys in ${avail_sys[@]}; do
            grep -Eq "^ $sys[0-9][0-9] " "$ctrl_file" || continue
            if ! grep -q "^+ .*$sys[0-9][0-9]" "$sp3"; then
                echo -e "$MSGWAR no satellite orbit for GNSS ($sys): $custom_pro_sp3"
            fi
        done
    fi

    # Check code/phase biases
    if [ -f "$fcb" ]; then
        if ! grep -q "^ OSB " "$fcb"; then
            echo -e "$MSGERR unsupported bias type (not OSB): $custom_pro_fcb"
            [ "$AR" == "Y" ] && return 1 || rm -f "$fcb"
        fi

        local obs_cand sys_num sys num
        local tna_cand tna fq1 fq2
        local freq_cmb=($(get_ctrl "$ctrl_file" "Frequency combination"))
        for sys_num in ${freq_cmb[@]}; do
            sys=${sys_num:0:1}
            grep -q "^ $sys[0-9][0-9] " "$ctrl_file" || continue
            obs_cand+=("${sys}XX")
            num=${sys_num:1:1}
            obs_cand+=("${sys}C${num}")
            obs_cand+=("${sys}L${num}")
            num=${sys_num:2:1}
            obs_cand+=("${sys}C${num}")
            obs_cand+=("${sys}L${num}")
        done

        for sys_num in ${obs_cand[@]}; do
            sys="${sys_num:0:1}"
            grep -q "OSB .* $sys[0-9][0-9] .* ${sys_num:1:2}" "$fcb" || continue
            obs_cand=("${obs_cand[@]/"$sys_num"}")
            obs_cand=("${obs_cand[@]/"${sys}XX"}")
        done

        local avail_sys=("G" "R" "E" "C" "J")
        for sys in ${avail_sys[@]}; do
            grep -Eq "^ $sys[0-9][0-9] " "$ctrl_file" || continue
            if echo "${obs_cand[@]}" | grep -q "${sys}XX"; then
                echo -e "$MSGWAR no OSB for GNSS ($sys)"
                continue
            fi
            sys_num=($(echo "${obs_cand[@]}" | grep -o "$sys\w\w" | cut -c 2-3))
            if [ ${#sys_num[@]} -ne 0 ]; then
                [[ "$sys" == "R" ]] && [[ "${sys_num[@]}" == "L1 L2" ]] && continue
                echo -e "$MSGWAR no OSB for GNSS ($sys): ${sys_num[@]}"
                continue
            fi
        done

        if  [[ $(head -1 "$rinexobs" | cut -c 6) == "3" ]]; then
            local obstypes=$(grep "SYS / # / OBS TYPES" "$rinexobs")
            while IFS= read -r line; do
                if [[ " " != $(echo "$line" | cut -c 1) ]]; then
                    sys=$(echo "$line" | cut -c 1) && tna_cand=("") && fq1=" " && fq2=" "
                    grep -q "^ $sys[0-9][0-9] " "$config"      || continue
                    echo "${obs_cand[@]}" | grep -q "${sys}XX" && continue
                    echo "${obs_cand[@]}" | grep -q "${sys}${fq1}X" || fq1=$(echo "${freq_cmb[@]}" | grep -o "${sys}[0-9][0-9]" | cut -c 2)
                    echo "${obs_cand[@]}" | grep -q "${sys}${fq2}X" || fq2=$(echo "${freq_cmb[@]}" | grep -o "${sys}[0-9][0-9]" | cut -c 3)
                fi
                tna_cand+=($(echo "$line" | grep -o " C${fq1}[A-Z] "))
                tna_cand+=($(echo "$line" | grep -o " C${fq2}[A-Z] "))
                for tna in ${tna_cand[@]}; do
                    grep -q "OSB .* $sys[0-9][0-9] .* $tna" "$fcb" && continue
                    echo -e "$MSGWAR no OSB for GNSS ($sys): $tna"
                    tna_cand=("${tna_cand[@]/$tna}")
                done
            done <<< "$obstypes"
        fi
    fi

    # APC model for code/phase biases
    local apc_setting=$(get_ctrl "$ctrl_file" "PCO on wide-lane")
    if [ "$apc_setting" == "Default" ]; then
        local apc_keyword=$(grep "APC_MODEL" "$fcb" | head -1 | awk '{print($2)}')
        if [[ -n "$apc_keyword" ]] && [[ "$apc_keyword" =~ ^NO* ]]; then
            sedi "/^PCO on wide-lane/s/ = .*/ = NO/"  "$ctrl_file"
        else
            sedi "/^PCO on wide-lane/s/ = .*/ = YES/" "$ctrl_file"
        fi
    fi

    # IGS ANTEX
    local abs_atx abs_url
    abs_atx="$(grep "SYS / PCVS APPLIED" $clk | head -1 | cut -c21-34 | tr 'A-Z' 'a-z' | sed 's/r3/R3/; s/ //g')"

    ### CODE MGEX ANTEX
    echo "$custom_pro_clk" | grep -qE "^ *(COD0MGX|COM)"
    if [[ $? -eq 0 ]] && [[ $abs_atx == "igs14" ]]; then
        [[ "$mjd_s" -le 59336 ]] && abs_atx="M14.ATX" || abs_atx="M20.ATX"
        atx_url="ftp://ftp.aiub.unibe.ch/CODE_MGEX/CODE/$abs_atx"
    fi

    if [ -n "$abs_atx" ]; then
        [[ "$abs_atx" =~ \.(ATX|atx)$ ]] || abs_atx="${abs_atx}.atx"
        echo -e "$MSGINF Prepare IGS ANTEX file: $abs_atx ..."
    else
        [[ "$OFFLINE" == "NO" ]] && abs_atx=$(curl https://files.igs.org/pub/station/general/ | grep -Eo "igs[0-9]{2}_[0-9]{4}.atx" | tail -1)
        [ -n "$abs_atx" ] || abs_atx=$(ls "$table_dir" | grep -Eo "igs[0-9]{2}_[0-9]{4}.atx" | tail -1)
        echo -e "$MSGINF Prepare IGS ANTEX file: $abs_atx ..."
        echo -e "$MSGWAR no PCO/PCV model specified in clock product $clk, use $table_dir/$abs_atx instead"
    fi

    if [ -f "$table_dir/$abs_atx" ]; then
        ln -sf "$table_dir/$abs_atx" abs_igs.atx
    else
        if [ -n "$atx_url" ]; then
            WgetDownload "$atx_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to download ANTEX file: $abs_atx"
                echo -e "$MSGINF please download from $atx_url to $table_dir for processing"
                return 1
            fi
        elif [[ $abs_atx =~ ^igs[0-9]{2} ]]; then
            atx_url="https://files.igs.org/pub/station/general/$abs_atx"
            WgetDownload "$atx_url"
            if [ $? -ne 0 ]; then
                atx_url="https://files.igs.org/pub/station/general/pcv_archive/$abs_atx"
                WgetDownload "$atx_url"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: failed to download ANTEX file: $abs_atx"
                    echo -e "$MSGINF please download from $atx_url to $table_dir for processing"
                    return 1
                fi
            fi
        elif [[ $abs_atx =~ ^igsR3 ]]; then
            atx_url="ftp://igs-rf.ign.fr/pub/IGSR3/$abs_atx"
            WgetDownload "$atx_url"
            if [ $? -ne 0 ]; then
                atx_url="ftp.aiub.unibe.ch/users/villiger/$abs_atx"
                WgetDownload "$atx_url"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: failed to download ANTEX file: $abs_atx"
                    echo -e "$MSGINF please download from $atx_url to $table_dir for processing"
                    return 1
                fi
            fi
        fi
        if [ -f "$abs_atx" ]; then
            [ -f "$table_dir/$abs_atx" ] || cp -f "$abs_atx" "$table_dir/"
            mv -f "$abs_atx" abs_igs.atx
        else
            echo -e "$MSGERR PrepareProducts: no IGS ANTEX file: $table_dir/$abs_atx"
            return 1
        fi
    fi

    ### LEO ANTEX 
    if [ "$mode" == "L" ]; then
        grep -q "^${site}  .*TYPE / SERIAL NO$" abs_igs.atx || leoatx.py abs_igs.atx ${site}
    fi

    echo -e "$MSGINF Prepare IGS ANTEX file: $abs_atx done"

    # Position SINEX solution 
    if [ "$mode" == "F" ]; then
        [ "$OFFLINE" == "NO" ] && mkdir -p "$product_ssc_dir"
        local wkdow=($(mjd2wkdow $mjd_s))
        if [ ${wkdow[0]} -lt 2238 ]; then
            local ssc="igs${ymd_s:2:2}P${wkdow[0]}${wkdow[1]}.ssc"
            local ssc_cmp="${ssc}.Z"
        else
            local ssc="IGS0OPSSNX_${ymd_s}${doy_s}0000_01D_01D_CRD.SNX"
            local ssc_cmp="${ssc}.gz"
        fi
        CopyOrDownloadProduct "$product_ssc_dir/$ssc"
        if [ $? -ne 0 ]; then
            for ssc_url in "ftp://igs.gnsswhu.cn/pub/gps/products/$wkdow" \
                           "ftp://nfs.kasi.re.kr/gps/products/$wkdow"     \
                           "ftp://gssc.esa.int/cddis/gnss/products/$wkdow"; do
                CopyOrDownloadProduct "$product_ssc_dir/$ssc_cmp" "$ssc_url/$ssc_cmp"
                [ $? -eq 0 ] && break
            done
            if [ ! -f "$ssc_cmp" ]; then
                echo -e "$MSGWAR PrepareProducts: failed to download position SINEX solution: $ssc_cmp"
                if [ ! -f sit.xyz ]; then
                    echo -e "$MSGERR no position SINEX solution: $ssc"
                    echo -e "$MSGINF please download from $ssc_url to $product_ssc_dir for processing"
                    return 1
                fi
            fi
        fi
        [ -f "$ssc_cmp" ] && gunzip -f "$ssc_cmp"
    fi

    # GIM (global ionospheric maps)
    local ion tec num
    if [ "$(get_ctrl "$config" "Iono 2nd")" == "YES" ]; then
        echo -e "$MSGSTA Downloading GIM products ..."
        [ "$OFFLINE" == "NO" ] && mkdir -p "$product_ion_dir"
        for mjd in $(seq $mjd_s $mjd_e); do
            local ydoy=($(mjd2ydoy $mjd))
            if [ $mjd -le 59909 ]; then
                local ion_tmp="CODG${ydoy[1]}0.${ydoy[0]:2:2}I"
                local ion_cmp="${ion_tmp}.Z"
            else
                local ion_tmp="COD0OPSFIN_${ydoy[0]}${ydoy[1]}0000_01D_01H_GIM.INX"
                local ion_cmp="${ion_tmp}.gz"
            fi
            local ion_url="ftp://ftp.aiub.unibe.ch/CODE/${ydoy[0]}/${ion_cmp}"
            CopyOrDownloadProduct "$product_ion_dir/$ion_tmp"
            if [ $? -ne 0 ]; then
                CopyOrDownloadProduct "$product_ion_dir/$ion_cmp" "$ion_url"
                if [ $? -ne 0 ]; then
                    echo -e "$MSGERR PrepareProducts: failed to download GIM product: $ion_cmp"
                    echo -e "$MSGINF please download from $ion_url to $product_ion_dir for processing"
                    return 1
                fi
            fi
            [ -f "$ion_cmp" ] && gunzip -f "$ion_cmp"
            ion="$ion $ion_tmp"
        done
        echo -e "$MSGSTA Downloading GIM products done"
        tec="tec_${ymd_s}${doy_s}"
        num="$(echo $ion | wc -w)"
        if [ $num -gt 1 ]; then
            MergeFiles "$(pwd)" "$ion" "$tec"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to merge GIM products: $ion -> $tec"
                return 1
            fi
            rm -f "$ion"
        else
            mv -f "$ion_tmp" "$tec"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to rename GIM product: $ion_tmp -> $tec"
                return 1
            fi
        fi
    fi

    # VMF (Vienna mapping function) grid
    local tmpy mjd hour

    grep '^ [0-9a-zA-Z]\{4\} .*VM1' "$config" &>/dev/null
    if [ $? -eq 0 ]; then
        echo -e "$MSGSTA Downloading VMF1 grid products ..."
        [ "$OFFLINE" == "NO" ] && mkdir -p "$product_vmf_dir"

        # Previous Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_s-1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        local vmf="VMFG_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H18"
        local vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$product_vmf_dir/$vmf"
        if [ $? -ne 0 ]; then
            CopyOrDownloadProduct "$product_vmf_dir/$vmf" "$vmf_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to download VMF1 grid product: $vmf"
                echo -e "$MSGINF please download from $vmf_url to $product_vmf_dir for processing"
                return 1
            fi
        fi

        # Current Day (for interpolation)
        for mjd in $(seq $mjd_s $mjd_e); do
            tmpy=($(mjd2ydoy $((mjd))))
            tmpy=($(ydoy2ymd ${tmpy[*]}))
            for hour in $(seq -w 00 06 18); do
                vmf="VMFG_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H${hour}"
                vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${tmpy[0]}/${vmf}"
                CopyOrDownloadProduct "$product_vmf_dir/$vmf"
                if [ $? -ne 0 ]; then
                    CopyOrDownloadProduct "$product_vmf_dir/$vmf" "$vmf_url"
                    if [ $? -ne 0 ]; then
                        echo -e "$MSGERR PrepareProducts: failed to download VMF1 grid product: $vmf"
                        echo -e "$MSGINF please download from $vmf_url to $product_vmf_dir for processing"
                        return 1
                    fi
                fi
            done
        done

        # Next Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_e+1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMFG_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H00"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$product_vmf_dir/$vmf"
        if [ $? -ne 0 ]; then
            CopyOrDownloadProduct "$product_vmf_dir/$vmf" "$vmf_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to download VMF1 grid product: $vmf"
                echo -e "$MSGINF please download from $vmf_url to $product_vmf_dir for processing"
                return 1
            fi
        fi

        cat VMFG_* > vmf_${ymd_s}${doy_s} || return 1
        echo -e "$MSGSTA Downloading VMF1 grid product done"
    fi

    grep '^ [0-9a-zA-Z]\{4\} .*VM3' "$config" &>/dev/null
    if [ $? -eq 0 ]; then
        echo -e "$MSGSTA Downloading VMF3 grid products ..."
        [ "$OFFLINE" == "NO" ] && mkdir -p "$product_vmf_dir"

        # Previous Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_s-1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        local vmf="VMF3_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H18"
        local vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$product_vmf_dir/$vmf"
        if [ $? -ne 0 ]; then
            CopyOrDownloadProduct "$product_vmf_dir/$vmf" "$vmf_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to download VMF3 grid product: $vmf"
                echo -e "$MSGINF please download from $vmf_url to $product_vmf_dir for processing"
                return 1
            fi
        fi

        # Current Day (for interpolation)
        for mjd in $(seq $mjd_s $mjd_e); do
            tmpy=($(mjd2ydoy $((mjd))))
            tmpy=($(ydoy2ymd ${tmpy[*]}))
            for hour in $(seq -w 00 06 18); do
                vmf="VMF3_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H${hour}"
                vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/${tmpy[0]}/${vmf}"
                CopyOrDownloadProduct "$product_vmf_dir/$vmf"
                if [ $? -ne 0 ]; then
                    CopyOrDownloadProduct "$product_vmf_dir/$vmf" "$vmf_url"
                    if [ $? -ne 0 ]; then
                        echo -e "$MSGERR PrepareProducts: failed to download VMF3 grid product: $vmf"
                        echo -e "$MSGINF please download from $vmf_url to $product_vmf_dir for processing"
                        return 1
                    fi
                fi
            done
        done

        # Next Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_e+1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMF3_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H00"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$product_vmf_dir/$vmf"
        if [ $? -ne 0 ]; then
            CopyOrDownloadProduct "$product_vmf_dir/$vmf" "$vmf_url"
            if [ $? -ne 0 ]; then
                echo -e "$MSGERR PrepareProducts: failed to download VMF3 grid product: $vmf"
                echo -e "$MSGINF please download from $vmf_url to $product_vmf_dir for processing"
                return 1
            fi
        fi

        cat VMF3_* > vmf_${ymd_s}${doy_s} || return 1
        echo -e "$MSGSTA Downloading VMF3 grid products done"
    fi

    # Rename products
    mv -f "${clk}" sck_${ymd_s}${doy_s} || return 1
    mv -f "${erp}" igserp               || return 1
    [ -e  "${fcb}" ] && mv -f "${fcb}" "fcb_${ymd_s}${doy_s}"
    [ -e  "${att}" ] && mv -f "${att}" "att_${ymd_s}${doy_s}"
    [ -e  "${ssc}" ] && mv -f "${ssc}" "snx_${ymd_s}${doy_s}"
    [ -e  "${lat}" ] && mv -f "${lat}" "lat_${ymd_s}${doy_s}"

    echo -e "$MSGSTA PrepareProducts done"
}

MergeFiles() { # purpose : merge multiple files into a single one
               # usage   : MergeFiles dir infile outfile
    local dir="$1"
    local infile="$2"
    local outfile="$3"
    rm -f "$outfile"
    for f in $infile; do
        cat "$dir/$f" >> "$outfile"
    done
    [ -f "$outfile" ] && return 0 || return 1
}

CopyOrDownloadProduct() { # purpose : copy or download a product
                          # usage   : CopyOrDownloadProduct copy url
    local copy="$1"
    local url="$2"
    local file=$(basename "$copy")

    # Try using cache from path "copy"
    if [ "$OFFLINE" = "YES" ] || [ "$USECACHE" = "YES" ]; then
        if [ -f "$copy" ]; then
            cp -f "$copy" .
        elif [ -f "$copy".gz ]; then
            cp -f "$copy".gz . && gunzip -f "$file".gz
        elif [ -f "$copy".Z ]; then
            cp -f "$copy".Z  . && gunzip -f "$file".Z
        fi
    fi

    # Try downloading file from "url"
    if [ ! -f "$file" ]; then
        WgetDownload "$url" || return 1
        if [ "$OFFLINE" = "YES" ] || [ "$USECACHE" = "YES" ]; then
            ls "$copy".@(gz|Z) &>/dev/null || cp -f "$(basename "$url")" "$copy"
        fi
    fi

    [ -f "$file" ] && return 0 || return 1
}

WgetDownload() { # purpose : download a file with wget
                 # usage   : WgetDownload url
    local url="$1"
    local arg="-q -nv -nc -c -t 3 --connect-timeout=10 --read-timeout=60"
    [ -n "$url" ] && [ "$OFFLINE" = "NO" ] || return 1

    wget --help | grep -q "\--show-progress" && arg="$arg --show-progress"
    local cmd="wget $arg $url"
    echo "$cmd" | bash

    [ -e $(basename "$url") ] && return 0  || return 1
}

LastYearMonth() { # purpose : get last year-month
                  # usage   : LastYearMonth year month
    local year=$1
    local mon=$((10#$2))
    [ $((mon-1)) -lt 1  ] && mon=12 && year=$((year-1)) || mon=$((mon-1))
    printf "%4d %02d\n" $year $mon
}

CleanAll() { # purpose : clean all files generated by PRIDE-PPPAR in the work directory
             # usage   : CleanAll year doy
    local year=$1
    local doy=$2
    rm -f igserp config\.*
    local types typ
    types=(rck ztd htg amb res stt cst neq att fcb orb sck snx)
    for typ in ${types[*]}; do
        rm -f ${typ}_${year}${doy}
    done
    types=(log pos kin)
    for typ in ${types[*]}; do
        rm -f ${typ}_${year}${doy}_${site}
    done
}

CleanMid() { # purpose : clean the intermediate files generated by PRIDE-PPPAR in the work directory
             # usage   : CleanMid year doy
    local year=$1
    local doy=$2
    local types typ
    types=(rck ztd htg amb res stt cst neq)
    for typ in ${types[*]}; do
        rm -f ${typ}_${year}${doy}
    done
    types=(log pos kin)
    for typ in ${types[*]}; do
        rm -f ${typ}_${year}${doy}_${site}
    done
}

Execute() {
    local cmd="$1"
    if [ $# -gt 1 ]; then
        local outp="$2"
    fi
    time=$(date +'%Y-%m-%d %H:%M:%S')
    if [ $# -gt 1 ]; then
        echo "$cmd" | bash > "$outp"
    else
        echo "$cmd" | bash
    fi
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} execution ok"
        return 0
    else
        echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} execution failed"
        return 1
    fi
}

ExecuteWithoutOutput() {
    local cmd="$1"
    time=$(date +'%Y-%m-%d %H:%M:%S')
    echo "$cmd" | bash &>/dev/null
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} execution ok"
        return 0
    else
        echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} execution failed"
        echo -e "$MSGINF Here is the output:\n"
        echo "$cmd" | bash
        return 1
    fi
}

snx2sit() {
    local site="$1"
    local mjd="$2"
    local wkdow=($(mjd2wkdow $mjd))
    local ydoy=($(mjd2ydoy $mjd))
    awk -v site=$(echo "$site" | tr 'a-z' 'A-Z') '{ 
            if ($1 == "+SOLUTION/ESTIMATE") fg = 1
            if ($1 == "-SOLUTION/ESTIMATE") fg = 0
            if (fg == 1) {
                if ($2 == "STAX") {
                   x = $9; sigx = $10;
                }
                if ($2 == "STAY") {
                   y = $9; sigy = $10;
                }
                if ($2 == "STAZ" && $3 == site) {
                    z = $9; sigz = $10;
                    printf(" %25.6f %25.6f %25.6f %25.6f %25.6f %25.6f\n", x, y, z, sigx, sigy, sigz)
                }
            } else {
                x=0.0; y=0.0; z=0.0; sigx=0.0; sigy=0.0; sigz=0.0;
            }
        }' "snx_${ydoy[0]}${ydoy[1]}" 
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

xyz2blh() {
    local x=$1
    local y=$2
    local z=$3
    echo "$x $y $z" | awk '
        BEGIN {
            F = 298.257223563;
            A = 6378137.0;
            B = A - A/F;
            E = 1 - (B/A)^2;
        } {
            x=$1; y=$2; z=$3;
            d = sqrt(x^2+y^2);
            h0 = sqrt(d^2+z^2) - A;
            b0 = z/d/(1-E*A/(A+h));
            while (i++ < 5) {
                n = A/sqrt(1-E*(1/(1+1/b0^2)));
                h = d/(1/sqrt(1+b0^2)) - n;
                b = z/d + n/(n+h) * E * b0;
                h0 = h;
                b0 = b;
            }
        } END {
            b = atan2(b, 1) * 180/atan2(0, -1);
            l = atan2(y, x) * 180/atan2(0, -1);
            if (l < 0) l += 360;
            printf("%15.7f%15.7f%15.4f\n", b, l, h)
        }'
}

######################################################################
##                               Entry                              ##
######################################################################

main "$@"
