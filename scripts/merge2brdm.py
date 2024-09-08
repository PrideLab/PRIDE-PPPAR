#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Merge GPS and GLONASS broadcast ephemerides to brdm file        ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 3.0                                                         ##
##                                                                           ##
##  DATE   : Mar-07, 2022                                                    ##
##                                                                           ##
##              @ GNSS RESEARCH CENTER, WUHAN UNIVERSITY, 2023               ##
##                                                                           ##
##    Copyright (C) 2022 by Wuhan University                                 ##
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

import fnmatch
import os
import sys

if len(sys.argv) != 3:
    print("usage: merge2brdm.py brdcddd.yyn brdcddd.yyg")
    sys.exit(1)

brdn = sys.argv[1]
brdg = sys.argv[2]

DDD = brdn[4:7]
YY  = brdn[9:11]
if brdg[4:7] == DDD and brdg[9:11] == YY:
    brdm = "brdm"+DDD+"0."+YY+"p"
else:
    print("***ERROR: inconsistent file name:")
    print("  "+brdn+" "+brdg)
    sys.exit(1)

fm = open(brdm,"w")
fm.write("     3.04           NAVIGATION DATA     M (Mixed)           RINEX VERSION / TYPE\n")

try:
    fn = open(brdn)
    lines = fn.readlines()
    inHeader = True
except:
    print("***ERROR: unable to open or read file "+brdn)
    sys.exit(1)

i = 1
while (i <= len(lines)):
    try:
        if (not inHeader):
            line = lines[i].replace("D","e")
            prn  = int(line[ 0: 2])
            yyyy = int(line[ 3: 5]) + 2000
            mm   = int(line[ 6: 8])
            dd   = int(line[ 9:11])
            hh   = int(line[12:14])
            mi   = int(line[15:17])
            ss   = round(float(line[18:22]))
            num2 = eval(line[22:41])
            num3 = eval(line[41:60])
            num4 = eval(line[60:79])
            fm.write("G{:02d} {:04d} {:02d} {:02d} {:02d} {:02d} {:02d}{: .12e}{: .12e}{: .12e}\n".format(
                prn, yyyy, mm, dd, hh, mi, int(ss), num2, num3, num4))
            for t in range(1,7):
                line = lines[i+t].replace("D","e")
                num1 = eval(line[ 3:22])
                num2 = eval(line[22:41])
                num3 = eval(line[41:60])
                num4 = eval(line[60:79])
                fm.write("    {: .12e}{: .12e}{: .12e}{: .12e}\n".format(num1, num2, num3, num4))
            line = lines[i+7].replace("D","e")
            num1 = eval(line[ 3:22])
            num2 = eval(line[22:41])
            fm.write("    {: .12e}{: .12e}\n".format(num1, num2))
            i = i + 8
            if (i >= len(lines)):
                break
        else:
            if ("PGM / RUN BY / DATE" == lines[i][60:79]):
                fm.write(lines[i])
            if ("LEAP SECONDS"        == lines[i][60:72]):
                leap_n = int(lines[i][1:6])
                fm.write(lines[i])
            if ("END OF HEADER"       == lines[i][60:73]):
                inHeader = False
                fm.write(lines[i])
            i = i + 1
    except:
        print("***ERROR: unexpected ERROR occured at line %d of file %s:" % (i, brdn))
        print(lines[i])
        sys.exit(1)

fn.close()

try:
    fg = open(brdg)
    lines = fg.readlines()
    inHeader = True
except:
    print("***ERROR: unable to open or read file "+brdg)
    sys.exit(1)

i = 1
while (i <= len(lines)):
    try:
        if (not inHeader):
            line = lines[i].replace("D","e")
            prn  = int(line[ 0: 2])
            yyyy = int(line[ 3: 5]) + 2000
            mm   = int(line[ 6: 8])
            dd   = int(line[ 9:11])
            hh   = int(line[12:14])
            mi   = int(line[15:17])
            ss   = round(float(line[18:22]))
            num2 = eval(line[22:41])
            num3 = eval(line[41:60])
            num4 = eval(line[60:79])
            fm.write("R{:02d} {:04d} {:02d} {:02d} {:02d} {:02d} {:02d}{: .12e}{: .12e}{: .12e}\n".format(
                prn, yyyy, mm, dd, hh, mi, int(ss), num2, num3, num4))
            for t in range(1,4):
                line = lines[i+t].replace("D","e")
                num1 = eval(line[ 3:22])
                num2 = eval(line[22:41])
                num3 = eval(line[41:60])
                num4 = eval(line[60:79])
                fm.write("    {: .12e}{: .12e}{: .12e}{: .12e}\n".format(num1, num2, num3, num4))
            i = i + 4
            if (i >= len(lines)):
                break
        else:
            if ("LEAP SECONDS"  == lines[i][60:72]):
                leap_g = int(lines[i][1:6])
            if ("END OF HEADER" == lines[i][60:73]):
                inHeader = False
            i = i + 1
    except:
        print("***ERROR: unexpected ERROR occured at line %d of file %s:" % (i, brdg))
        print(lines[i])
        sys.exit(1)

fg.close()
fm.close()
