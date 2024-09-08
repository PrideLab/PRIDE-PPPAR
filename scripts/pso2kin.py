#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Convert PSO file to kin format of PRIDE PPP-AR                  ##
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

import os
import sys
from datetime import datetime
from datetime import timedelta
import platform

global PATHSEP
if platform.system() == "Windows":
    PATHSEP = "\\"
else:
    PATHSEP = '/'


# time functions
def mjd2ydoy(mjd):
    time = mjd2time(mjd)
    time0 = datetime(time.year, 1, 1, 0, 0, 0)
    return time.year, (time - time0).days + 1


def mjd2time(mjd):
    time0 = datetime(1858, 11, 17, 0, 0, 0, 0)
    time1 = time0 + timedelta(days=mjd)
    return time1


def secG_ymd(sec):
    time0 = datetime(2000, 1, 1, 12, 00, 00)
    time1 = time0 + timedelta(seconds=sec)
    return time1


def time2mjd(dateT):
    t0 = datetime(1858, 11, 17, 0, 0, 0, 0)
    mjd = (dateT-t0).days
    sod = dateT.hour*3600.0+dateT.minute*60.0+dateT.second+dateT.microsecond/1000000.0
    return mjd, sod


def t2mjd(t):
    time0 = datetime(1858, 11, 17, 0, 0, 0, 0)
    time = datetime(t[0], t[1], t[2], 0, 0, 0)
    mjd = (time-time0).days
    sod = t[3]*3600+t[4]*60+t[5]
    return mjd, sod


# find & read pso files
def find_allpso(leodir, leoid, mjds, mjde):
    psofiles = []

    for mjd in range(mjds, mjde+1):
        year, doy = mjd2ydoy(mjd)
        psofl = leodir + PATHSEP + "pso_%4d%03d_%s" % (year, doy, leoid)

        if not os.path.exists(psofl):
            print("pso2kin.py Warning: No such pso file: %s!\n" % psofl)
            sys.exit(1)
        if not os.access(psofl, os.R_OK):
            print("pso2kin.py Warning: Pso file %s is unreadable!\n" % psofl)
            sys.exit(1)

        psofiles.append(psofl)

    return psofiles


def readpso(psofiles, leoid):
    if 'gra' in leoid:
        data = read_grapso(psofiles)
    else:
        print("pso2kin.py Warning: unsupported LEO satellite!")
        sys.exit(1)
    return data


def read_grapso(psofiles):
    data = []
    for psofl in psofiles:
        try:
            with open(psofl, 'r') as rdpso:
                while (1):
                    line = rdpso.readline().strip("\n")
                    if "# End of YAML header" in line or "END OF HEADER" in line:
                        break
                lines = rdpso.readlines()

            for line in lines:
                strs = line.split()
                t = secG_ymd(int(strs[0]))
                mjd, sod = time2mjd(t)
                x = float(strs[3])
                y = float(strs[4])
                z = float(strs[5])
                data.append([mjd, sod, x, y, z])
        except OSError:
            print("pso2kin.py Warning: open pso file %s failed!" % psofl)
            sys.exit(1)
    return data


# convert pso files to kin_ format
def pso2kin(data, leoid, mjds):
    year0, doy0 = mjd2ydoy(mjds)
    kinfl = "kin_%4d%03d_%s" % (year0, doy0, leoid)

    try:
        with open(kinfl, 'w') as wrkin:
            wrkin.write("Kinematic Trajectory          GRAC                          COMMENT\n")
            wrkin.write("                                                            INTERVAL\n")
            wrkin.write("                                                            END OF HEADER\n")
            for i in range(len(data)):
                mjd = int(data[i][0])
                [sod, x, y, z] = [float(data[i][j+1]) for j in range(4)]
                wrkin.write("{}{:> 10.2f}{:> 15.4f}{:> 14.4f}{:> 14.4f}\n" .format(mjd, sod, x, y, z))
    except IOError:
        print("pso2kin.py Warning : can not open kin_ file %s!" % kinfl)
        sys.exit(1)
    '''
    mjd1 = int(data[-1][0])
    sod1 = float(data[-1][1])
    [x0, y0, z0] = [float(data[0][j+2]) for j in range(3)]
    print("Position :%17.4f%17.4f%17.4f%17.4f\n" % (x, y, z))
    print("Duration : %4d %02d %02d %02d %02d  %4.2f%11.2f" % (year0)) 
    # Should be determined by the RINEXo file (spp?)
    '''


def main(leodir, leoid, mjds, mjde):
    psofiles = find_allpso(leodir, leoid, mjds, mjde)
    data = readpso(psofiles, leoid)
    pso2kin(data, leoid, mjds)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage  : pso2kin.py leodir leoid mjds mjde") 
        sys.exit(1)
    leodir = sys.argv[1]
    leoid = sys.argv[2]
    mjds = int(sys.argv[3])
    mjde = int(sys.argv[4])
    main(leodir, leoid, mjds, mjde)
