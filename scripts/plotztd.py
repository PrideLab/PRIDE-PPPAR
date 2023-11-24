#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Plot zenith troposphere delays                                  ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 2.2                                                         ##
##                                                                           ##
##  DATE   : Jan-18, 2022                                                    ##
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

import sys
import os
import math
import warnings
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
 
def read_ztdfil(ztd_file):
    ztd = []
    with open(ztd_file, 'r') as f:
        ## read header
        site = "NONE"
        while (1):
            line = f.readline().strip("\n")
            if "STATION" in line:
                site = line[0:4]
            elif "END OF HEADER" in line:
                break
        ## output file
        outfile = "ztd_tmp"
        wr=open(outfile, 'w')
        if line == "":
            f.close()
        else:
            while (1):
                line = f.readline().strip("\n")
                if line == "":
                    break
                elif line[0] == "*":
                    continue
                else:
                    wr.write(site+" "+line+"\n")
        f.close
    wr.close

## main
if len(sys.argv) != 3:
    print('#usage  : plotztd.py ztd_filename png_filename')
    print('#example: plotztd.py ztd_2021149_rov1 ztd_2021149_rov1')
    sys.exit(0)

ztdflname = sys.argv[1]

if (os.path.isfile('ztd_tmp') == False):
    if (os.path.isfile(ztdflname) == True):
        read_ztdfil(ztdflname)
    else:
        print('No such file:' + ztdflname)
        sys.exit(1)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    data = np.loadtxt('ztd_tmp', dtype=str, comments='#')

if (len(data) == 0):
    print('Empty input file: ztd_tmp')
    sys.exit(1)

n = len(data)
dates_tmp = [0]*n
site = data[0,0]
for i in range(n):
    if data[i,6] == "60.000000":
        data[i,5] = int(data[i,5])+1
        data[i,6] = "0.000000"
    if data[i,5] == "60" :
        data[i,4] = int(data[i,4])+1
        data[i,5] = "0"
    dates_tmp[i] = datetime.datetime(int(data[i,1]), int(data[i,2]), int(data[i,3]),
            int(data[i,4]), int(data[i,5]),0) + datetime.timedelta(seconds=float(data[i,6]))

dts = np.zeros(n-1)
for i in range(n-1):
    dts[i] = (dates_tmp[i+1] - dates_tmp[i]).total_seconds()
intv = np.nanmin(dts)
dt = (dates_tmp[n-1]-dates_tmp[0]).total_seconds()
nmax = int(math.ceil(dt/intv)+1)
dates = [0]*nmax
ztds = np.zeros(nmax)
for i in range(n):
    dt = (dates_tmp[i] - dates_tmp[0]).total_seconds()
    if dt%intv < 1E-4 or abs(dt%intv - intv) < 1E-4:
        iepo = int(round(dt/intv))
        ztds[iepo] = float(data[i,7]) + float(data[i,8]) + float(data[i,9])

for i in range(nmax):
    dt = i*intv
    dates[i] = dates_tmp[0] + datetime.timedelta(seconds=dt)
    if ztds[i] == 0.0:
        ztds[i] = np.nan

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12,4))
fig.subplots_adjust(hspace=0.05)
fig.align_ylabels()

ax.plot(dates[:], ztds[:], linewidth=1.5, color='k')
ax.text(0.98,0.94, site, horizontalalignment='right', verticalalignment='top', 
fontsize=16, transform=ax.transAxes)

fwidth = 1.5
ax.spines['bottom'].set_linewidth(fwidth)
ax.spines['left'].set_linewidth(fwidth)
ax.spines['top'].set_linewidth(fwidth)
ax.spines['right'].set_linewidth(fwidth)

beg = dates[0]
end = dates[nmax-1]
dt = (end-beg).total_seconds()
nday = (end-beg).days+1

if dt < 86400:
    mint_intv = int(15*(math.ceil(dt/3600/3)))
    intvs = mdates.MinuteLocator(interval=mint_intv)
    if mint_intv >=60:
       hour_intv=int(math.ceil((mint_intv)/60))
       intvs = mdates.HourLocator(interval=hour_intv)
else:
    nintv = int(12/nday)
    hour_intv = int(math.ceil(24/nintv))
    if hour_intv%2 != 0:
        hour_intv = hour_intv+1
    intvs = mdates.HourLocator(interval=hour_intv)

ax.set_xlim(beg,end)
ax.tick_params(axis='both', which='both', labelsize=16, pad=4)
ax.tick_params(axis='both', which='major', length=6)
ax.tick_params(axis='both', which='minor', length=4)
ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
ax.xaxis.set_major_locator(intvs)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)

## multiple days
if nday > 1:
    days = mdates.DayLocator(interval=1)
    newax = ax.twiny()
    newax.set_xlim(beg,end)
    newax.set_frame_on(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.tick_params(axis='both', which='major', direction='in', length=8)
    newax.tick_params(axis='both', which='minor', direction='in', length=0)
    newax.tick_params(axis='both', which='both', labelsize=16, pad=20)
    newax.xaxis.set_major_locator(days)
    newax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

stdv = np.nanstd(ztds[:])
min_elev = np.nanmin(ztds[:])-2*stdv
max_elev = np.nanmax(ztds[:])+2*stdv
ax.set_ylim([min_elev, max_elev])

xlabel = "Time"
ylabel = "Zenith troposphere delays (m)"
ax.set_ylabel(ylabel, labelpad=3)
if nday > 1:
    ax.set_xlabel(xlabel, labelpad=23)
else:
    ax.set_xlabel(xlabel, labelpad=3)

pngflname = sys.argv[2]+'.png'
plt.savefig(pngflname, pad_inches=None, bbox_inches='tight')
plt.show()
os.system('rm ztd_tmp')
sys.exit(0)
