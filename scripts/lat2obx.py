#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Convert LEO quaternions file to ORBEX format                    ##
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

import sys
from datetime import datetime, timedelta

# time = [y, m, d, h, min, sec]
def ymd2doy(time):  
    time0 = datetime(time[0], 1, 1, 0, 0, 0)
    time1 = datetime(time[0], time[1], time[2], 0, 0, 0)
    return (time1-time0).days+1


# GRACE/GRACE-FO time stamp conversion
def gfot_time(grasec):
    time0 = datetime(2000, 1, 1, 12, 00, 00)
    time1 = time0 + timedelta(seconds=grasec)
    fsec = grasec-int(grasec)
    sec = time1.second+fsec
    time = [time1.year, time1.month, time1.day, time1.hour, time1.minute, sec]
    return time


def read_attfl(attfl, leocode):
    with open(attfl, 'r') as rdatt:
        while 1:
            line = rdatt.readline().strip('\n')
            if '# End of YAML header' in line:
                break
            if 'END OF HEADER' in line:
                break
            if 'Record' in line:
                break
        datatmp = rdatt.readlines()
    y, doy = get_ydoy(datatmp[0], leocode)
    return datatmp, y, doy


def get_ydoy(line, leocode):
    strs = line.split()
    if leocode[0:3] == 'gra':
        qtime = gfot_time(int(strs[0]))
    doy = ymd2doy(qtime)
    return qtime[0], doy


def line2q(line, leocode):
    strs = line.split()
    if leocode[0:3] == 'gra':
        qtime = gfot_time(int(strs[0]))
        q = [float(strs[i+3]) for i in range(4)]
    return qtime, q


def wr_obxh(wr, leocode):
    ref_fra = 'ECEF'
    if leocode[0:3] == 'gra':
        ref_fra = 'ECI'
    wr.write("%=ORBEX  0.09\n")
    wr.write("%%\n")
    wr.write("+FILE/DESCRIPTION\n")
    wr.write("%s%s\n" % (" DESCRIPTION".ljust(21, ' '), "LEO satellite attitude quaternions "))
    wr.write("%s%s\n" % (" CONVERTED_BY".ljust(21, ' '), "PRIDE Lab, GNSS Research Center, Wuhan University"))
    wr.write("%s%s\n" % (" CONTACT".ljust(21, ' '), "pride@whu.edu.cn"))
    wr.write("%s%s\n" % (" LEO_SATELLITE".ljust(21, ' '), leocode))
    wr.write("%s%s\n" % (" FRAME_TYPE".ljust(21, ' '), ref_fra))
    wr.write("-FILE/DESCRIPTION\n")
    wr.write("+EPHEMERIS/DATA\n")


def wr_obx(obxfl, leocode, datatmp):
    with open(obxfl, 'w') as wrobx:
        wr_obxh(wrobx, leocode)
        nepo = len(datatmp)
        for iepo in range(nepo):
            epoch, quat = line2q(datatmp[iepo], leocode)
            wrobx.write("%s%04d %02d %02d %02d %02d %15.12f %d\n" % ('## ', epoch[0], epoch[1],
                        epoch[2], epoch[3], epoch[4], epoch[5], 1))
            wrobx.write("{}{} {:> 17.16f} {:> 17.16f} {:> 17.16f} {:> 17.16f}\n" .format(" ATT".ljust(22, ' '), 4,
                        quat[0], quat[1], quat[2], quat[3]))
        wrobx.write("-EPHEMERIS/DATA\n")


def main(attfl, leocode):
    datatmp, y, day = read_attfl(attfl, leocode)
    doy = "%03d" % day
    obxfl = 'lat_' + str(y) + str(doy) + '_' + leocode
    wr_obx(obxfl, leocode, datatmp)


if len(sys.argv) != 3:
    print("Usage:	lat2obx.py latfl leocode ")
    print(" latfl       LEO attitude file")
    print(" leocode     LEO satellite name, format:xxxx")
    print("             GRACE/GRACE-FO :  graa/grab/grac/grad")
    sys.exit(0)


latfl = sys.argv[1]
leocode = sys.argv[2]
main(latfl, leocode)
