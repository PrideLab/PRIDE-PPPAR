#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Plot trajectory                                                 ##
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
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

def read_kinfil(kin_file):
    # read file
    site = "NONE"
    with open(kin_file, 'r') as f :
        ## header
        while (1):
            line = f.readline().strip("\n")
            if "STATION" in line:
                site = line[0:4]
            if "END OF HEADER" in line:
                break
        ## output file
        outfile = "kin_tmp"
        wr = open(outfile, 'w')
        if line == "":
            f.close()
        else:
            while (1):
                line = f.readline().strip("\n")
                if line == "":
                    break
                elif "*" in line:
                    continue
                else:
                    wr.write(site+" "+line+"\n")
            f.close()
    wr.close()

def mjd2time(mjd):
   t0 = datetime.datetime(1858,11,17,0,0,0,0)
   return t0+datetime.timedelta(days=mjd)

if len(sys.argv) != 3:
  print('#usage  : plottrack.py kin_filename png_filename')
  print('#example: plottrack.py kin_2021149_rov1 rov1_2021149')
  sys.exit(0)

kinflname = sys.argv[1]
pngflname = sys.argv[2]

if (os.path.isfile('kin_tmp') == False):
    if (os.path.isfile(kinflname) == True):
        read_kinfil(kinflname)
    else:
        print('Error: no such file: '+kinflname)
        sys.exit(1)
else:
    print('Warning: use old kin_tmp file')

## Prepare Date
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    datatmp = np.loadtxt('kin_tmp', dtype=str, comments='#')

n = len(datatmp)
if n==0: 
    print('Error: empty input file: kin_tmp')
    sys.exit(0)

site = datatmp[0,0]
data = np.zeros((n,5)) # mjd sod b l h
for i in range(n):
    data[i,0] = int(datatmp[i,1])
    data[i,1] = float(datatmp[i,2])
    for j in range(3):
        data[i,j+2] = float(datatmp[i,j+6])

dts = np.zeros(n-1)
for i in range(n-1):
    dts[i] = (data[i+1,0]-data[i,0])*86400.0+(data[i+1,1]-data[i,1])
intv = np.nanmin(dts)
dt = (data[n-1,0]-data[0,0])*86400.0+(data[n-1,1]-data[0,1])
nmax = int(math.ceil(dt/intv)+1)
blh = np.zeros((nmax,3))
dates = [0]*nmax
for i in range(n):
    dt = (data[i,0]-data[0,0])*86400.0+(data[i,1]-data[0,1])
    if dt%intv < 1E-4 or abs(dt%intv - intv) < 1E-4:
        iepo = int(round(dt/intv))
        for j in range(3):
            blh[iepo,j] = data[i,j+2]

beg = mjd2time(data[0,0])+datetime.timedelta(seconds=data[0,1])
for i in range(nmax):
    dt = i*intv
    time = beg+datetime.timedelta(seconds=dt)
    dates[i] = time
    if blh[i,0] == 0.0 and blh[i,1] == 0.0 and blh[i,2] == 0.0:
        for j in range(3):
            blh[i,j] = np.nan

## Plot
fig = plt.figure(figsize=(12,10))
gs = GridSpec(3, 1, figure=fig)
axs = []
axs.append(fig.add_subplot(gs[0:2]))
axs.append(fig.add_subplot(gs[2]))
fig.align_ylabels()

colors = ('k', 'r', 'b')
## plot blh
axs[0].plot(blh[:,1], blh[:,0], linewidth=1.5, color=colors[0])
axs[0].text(0.98, 0.97, site, horizontalalignment='right', verticalalignment='top', 
fontsize=18, transform=axs[0].transAxes)
axs[1].plot(dates[:], blh[:,2], linewidth=1.5, color=colors[0])

fwidth = 1.5
for ax in axs:
    ax.spines['bottom'].set_linewidth(fwidth)
    ax.spines['left'].set_linewidth(fwidth)
    ax.spines['top'].set_linewidth(fwidth)
    ax.spines['right'].set_linewidth(fwidth)

beg = mjd2time(data[0,0])+datetime.timedelta(seconds=data[0,1])
end = mjd2time(data[n-1,0])+datetime.timedelta(seconds=data[n-1,1])
dt = (end-beg).total_seconds()
mjd = np.unique(data[:,0])

if dt < 86400:
    mint_intv = int(15*(math.ceil(dt/3600/3)))
    intvs = mdates.MinuteLocator(interval=mint_intv)
    if mint_intv >=60:
       hour_intv=int(math.ceil((mint_intv)/60))
       intvs = mdates.HourLocator(interval=hour_intv)
else:
    nintv = int(12/len(mjd))
    hour_intv = int(math.ceil(24/nintv))
    if hour_intv%2 != 0:
        hour_intv = hour_intv+1
    intvs = mdates.HourLocator(interval=hour_intv)

axs[1].set_xlim(beg,end)
for ax in axs:
    ax.tick_params(axis='both', which='both', labelsize=18, pad=4)
    ax.tick_params(axis='both', which='major', length=6)
    ax.tick_params(axis='both', which='minor', length=4)
    ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
    ax.xaxis.label.set_size(18)
    ax.yaxis.label.set_size(18)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

axs[0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axs[1].xaxis.set_major_locator(intvs)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

## multiple days
if len(mjd) > 1:
    days = mdates.DayLocator()
    newax = axs[1].twiny()
    newax.set_xlim(beg,end)
    newax.set_frame_on(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.tick_params(axis='x', which='major', direction='in', length=8)
    newax.tick_params(axis='x', which='minor', direction='in', length=0)
    newax.tick_params(axis='both', which='both', labelsize=18, pad=20)
    newax.xaxis.set_major_locator(days)
    newax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

maxv = np.nanmax(blh[:,1])+np.nanstd(blh[:,1])
minv = np.nanmin(blh[:,1])-np.nanstd(blh[:,1])
axs[0].set_xlim([minv, maxv])

maxv = np.nanmax(blh[:,0])+np.nanstd(blh[:,0])
minv = np.nanmin(blh[:,0])-np.nanstd(blh[:,0])
axs[0].set_ylim([minv, maxv])

maxv = np.nanmax(blh[:,2])+np.nanstd(blh[:,2])
minv = np.nanmin(blh[:,2])-np.nanstd(blh[:,2])
axs[1].set_ylim([minv, maxv])

xlabels = ("Longitude (degree)", "Time")
ylabels = ("Latitude (degree)", "Height (m)")

axs[0].set_xlabel(xlabels[0], labelpad=3)
if len(mjd) > 1:
    axs[1].set_xlabel(xlabels[1],labelpad=23)
else:
    axs[1].set_xlabel(xlabels[1],labelpad=3)

for i in range(2):
    axs[i].set_ylabel(ylabels[i],labelpad=3)

plt.tight_layout()
plt.savefig(pngflname, pad_inches=None, bbox_inches='tight')
plt.show()
os.system('rm kin_tmp')
sys.exit(0)
