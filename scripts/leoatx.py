#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Add LEO APC information to .atx file as ANTEX format            ##
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


def mjd2time(mjd):
    time0 = datetime(1858, 11, 17, 0, 0, 0, 0)
    time1 = time0+timedelta(days=mjd)
    return time1


class ATX:
    def __init__(self):
        self.PCO = {'graa': {(52368, 58058): {'G01': "-0.4 -0.4 -451.142", 'G02': "-0.4 -0.4 -475.65"}},
                    'grab': {(52368, 58058): {'G01': "0.6 0.754 -451.73", 'G02': "0.6 0.754 -475.96"}},
                    'grac': {(58260, 99999): {'G01': "-1.28 260.23 -486.24", 'G02': "-1.28 260.23 -492.24"}},
                    'grad': {(58260, 99999): {'G01': "-1.07 260.04 -485.45", 'G02': "-1.07 260.04 -491.45"}},
                    }


def wr_leoatx(PCO, atx):
    try:
        PCO = ATX()
        with open(atx, 'a') as wr:
            for LEO in PCO.PCO.keys():
                for tspan in PCO.PCO[LEO].keys():
                    [mjds, mjde] = tspan
                    mpco = PCO.PCO[LEO][tspan]
                    stime = mjd2time(mjds)
                    ts = str(stime.year).rjust(6, ' ') + str(stime.month).rjust(6, ' ') + str(stime.day).rjust(6, ' ')
                    ts = ts + '0'.rjust(6, ' ') + '0'.rjust(6, ' ') + '0.0000000'.rjust(13, ' ')
                    lwrte = False
                    if mjde != 99999:
                        lwrte = True
                        etime = mjd2time(mjde)
                        te = str(etime.year).rjust(6, ' ') + str(etime.month).rjust(6, ' ') + str(etime.day).rjust(6, ' ')
                        te = te + '59'.rjust(6, ' ') + '59'.rjust(6, ' ') + '59.0000000'.rjust(13, ' ')
                    wr.write("\n%s\n" % "START OF ANTENNA".rjust(76, ' '))
                    wr.write("%s%s\n" % (LEO, "TYPE / SERIAL NO".rjust(76 - len(LEO))))
                    wr.write("%s%s\n" % ('     0.0'.ljust(60), 'DAZI'))
                    wr.write("%s%s\n" % ('     0.0  13.0   1.0'.ljust(60), 'ZEN1 / ZEN2 / DZEN'))
                    frenum = len(mpco)
                    wr.write("     %s%s\n" % (frenum, "# OF FREQUENCIES".rjust(76-5-len(str(frenum)))))
                    wr.write("%s%s\n" % (ts.ljust(60, ' '), 'VALID FROM'))
                    if lwrte:
                        wr.write("%s%s\n" % (te.ljust(60, ' '), 'VALID UNTIL'))
                    wr.write("%s%s\n" % (atx, "SINEX CODE".rjust(70 - len(atx))))
                    wr.write("%-60s%s\n" % ("LEO APC by leoatx.py", "COMMENT"))
                    for FRE in mpco:
                        pco = mpco[FRE]
                        wr.write("   %s%s\n" % (FRE, "START OF FREQUENCY".rjust(78-3-len(FRE))))
                        wr.write("    %s%s\n" % (pco, "NORTH / EAST / UP".rjust(77-4-len(pco))))
                        wr.write("%s%s\n" % ("   NOAZI", 14*"    0.00"))
                        wr.write("   %s%s\n" % (FRE, "END OF FREQUENCY".rjust(76-3-len(FRE))))
                    wr.write("%s\n" % "END OF ANTENNA".rjust(74, ' '))
    except OSError:
        print("  ###Warning(%s): Cannot write APC information of %s to %s!" % (sys.argv[0], LEO, atx))

    
def main(atx, LEO):
    PCO = ATX()
    if LEO not in PCO.PCO.keys():
        print("  ###Warning(%s):  No APC information for %s!" % (sys.argv[0], LEO))
        sys.exit(1)
    wr_leoatx(PCO, atx)


if __name__ == "__main__":
    if len(sys.argv) != 3 :
        print("Usage	: leoatx.py ANTEX LEOsat")
        print("Example	: leoatx.py abs_igs.atx grac")
        sys.exit(1)

    atx = sys.argv[1]
    LEO = sys.argv[2]
    main(atx, LEO)
