#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Plot kinematic positions                                        ##
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

def mjd2time(mjd):
    t0 = datetime.datetime(1858,11,17,0,0,0,0)
    return t0+datetime.timedelta(days=mjd)

if len(sys.argv) != 3 and len(sys.argv) != 6:
    print('#usage  : plotkin.py kin_filename png_filename [x_ref y_ref z_ref]')
    print('#example: plotkin.py kin_2021149_rov1 rov1_2021149')
    print('#if no x_ref y_ref z_ref, default x_avg y_avg z_avg')
    sys.exit(0)

kinflname = sys.argv[1]
pngflname = sys.argv[2]
if len(sys.argv) == 6:
    x_ref = sys.argv[3]
    y_ref = sys.argv[4]
    z_ref = sys.argv[5]

if len(sys.argv) == 3:
    os.system('xyz2enu ' + kinflname + ' enu_tmp')
else:
    os.system('xyz2enu ' + kinflname + ' enu_tmp' + ' ' + x_ref + ' ' + y_ref + ' ' + z_ref)

# Prepare Date
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    datatmp = np.loadtxt('enu_tmp', dtype=str, comments='#')

n = len(datatmp)
if n == 0: 
    print('Error: empty input file: enu_tmp')
    sys.exit(0)
    
data = np.zeros((n,7)) # mjd sod e n u PDOP Satnum
for i in range(n):
    data[i,0] = int(datatmp[i,0])
    for j in range(4):
        data[i,j+1] = float(datatmp[i,j+1])
    data[i,5] = int(datatmp[i,5])
    data[i,6] = float(datatmp[i,12])

dts = np.zeros(n-1)
for i in range(n-1):
    dts[i] = (data[i+1,0]-data[i,0])*86400.0+(data[i+1,1]-data[i,1])
intv = np.nanmin(dts)
dt = (data[n-1,0]-data[0,0])*86400.0+(data[n-1,1]-data[0,1])
nmax = int(math.ceil(dt/intv)+1)
enu = np.zeros((nmax,5))
dates = [0]*nmax
for i in range(n):
    dt = (data[i,0]-data[0,0])*86400.0+(data[i,1]-data[0,1])
    if dt%intv < 1E-4 or abs(dt%intv - intv) < 1E-4:
        iepo = int(round(dt/intv))
        for j in range(3):
            enu[iepo,j] = data[i,j+2]*100
        enu[iepo,3] = data[i,5]
        enu[iepo,4] = data[i,6]

beg = mjd2time(data[0,0])+datetime.timedelta(seconds=data[0,1])
for i in range(nmax):
    dt = i*intv
    time = beg+datetime.timedelta(seconds=dt)
    dates[i] = time
    if enu[i,0] == 0.0 and enu[i,1] == 0.0 and enu[i,2] == 0.0 \
            and enu[i,3] == 0.0 and enu[i,4] == 0.0:
        for j in range(5):
            enu[i,j] = np.nan

# Plot
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12,12))
fig.subplots_adjust(hspace=0.12)

colors = ('k', 'r', 'b')
strs = ('East ', 'North ', 'Up ')
# plot enu displacements
for i in range(3):
    axs[i].plot(dates[:], enu[:,i], linewidth=1.5, color=colors[0])
    rms = np.sqrt(np.nanmean(np.square(enu[:,i])))
    str1 = strs[i]+"%6.2f " %(rms)+'cm'
    axs[i].text(0.98, 0.92, str1, horizontalalignment='right', verticalalignment='top', 
                fontsize=16, transform=axs[i].transAxes)
# plot PDOP
lns1 = axs[3].plot(dates[:], enu[:,4], linewidth=2, color=colors[1], label='PDOP')
## plot Satellite number
ax2 = axs[3].twinx()
lns2 = ax2.plot(dates[:], enu[:,3], linewidth=2, color=colors[2], label='Satellite number')
lns = lns1+lns2
labels = [l.get_label() for l in lns]
axs[3].legend(lns, labels[:], ncol=2, columnspacing=1, edgecolor='black', fontsize=16)

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

axs[3].set_xlim(beg, end)
for ax in axs:
    ax.tick_params(axis='both', which='both', labelsize=16, pad=4)
    ax.tick_params(axis='both', which='major', length=6)
    ax.tick_params(axis='both', which='minor', length=4)
    ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
    ax.xaxis.set_major_locator(intvs)
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# multiple days
if len(mjd) > 1:
    days = mdates.DayLocator()
    newax = axs[3].twiny()
    newax.set_xlim(beg,end)
    newax.set_frame_on(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.tick_params(axis='x', which='major', direction='in', length=8)
    newax.tick_params(axis='x', which='minor', direction='in', length=0)
    newax.tick_params(axis='both', which='both', labelsize=16, pad=20)
    newax.xaxis.set_major_locator(days)
    newax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

maxv = math.ceil(np.nanmax(enu[:,0:2]))*1.5
minv = math.floor(np.nanmin(enu[:,0:2]))-1
for i in range(2):
    axs[i].set_ylim([minv, maxv])

maxv = math.ceil(np.nanmax(enu[:,2]))*1.5
minv = math.floor(np.nanmin(enu[:,2]))-1
axs[2].set_ylim([minv, maxv])

maxv = math.ceil(np.nanmax(enu[:,4]))*1.5
axs[3].set_ylim(0.5, maxv)

maxv = math.ceil(np.nanmax(enu[:,3]))+5
minv = math.floor(np.nanmin(enu[:,3]))-0.5
ax2.set_ylim(minv, maxv)

axs[3].set_xlim(beg,end)

axs[3].tick_params(axis='both', which='both', direction='in', right=False, top=True)
axs[3].tick_params(axis='y', which='both', colors='r')
ax2.tick_params(axis='both', labelsize=16, pad=4)
ax2.tick_params(axis='both', which='major', length=6)
ax2.tick_params(axis='both', which='minor', length=3)
ax2.tick_params(axis='both', which='both', direction='in', left=False, top=True)
ax2.tick_params(axis='y', which='both', colors='b')
ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))
ax2.yaxis.label.set_size(16)

xlabel = "Time"
ylabels = ('East (cm)', 'North (cm)', 'Up (cm)', "PDOP", "Satellite number")
if len(mjd) > 1:
    axs[3].set_xlabel(xlabel,labelpad=23)
else:
    axs[3].set_xlabel(xlabel,labelpad=3)
for i in range(3):
    axs[i].set_ylabel(ylabels[i], labelpad=3)
axs[3].set_ylabel(ylabels[3], labelpad=3, color='r')
axs[3].spines['left'].set_color('r')
axs[3].spines['right'].set_color('b')
ax2.spines['left'].set_color('r')
ax2.spines['right'].set_color('b')
ax2.set_ylabel(ylabels[4],labelpad=3, color='b')

plt.savefig(pngflname, pad_inches=None, bbox_inches='tight')
plt.show()
os.system('rm enu_tmp')
sys.exit(0)
