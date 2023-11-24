#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Plot pseudorange and carrier-phase residuals                    ##
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

def read_resfil(res_file):
    # read file
    sat = []
    wr = []
    with open(res_file, 'r') as f :
        ## header
        while (1):
            line = f.readline().strip("\n")
            if "END OF HEADER" in line:
                break
        ## output file
        outfile = "res_tmp"
        wr = open(outfile, 'w')
        if line == "":
            f.close()
        else:
            while (1):
                line = f.readline().strip("\n")
                if line == "":
                    break
                elif "CST" in line:
                    continue
                elif "TIM" in line:
                    time = line[32:]
                    continue
                odom = line.split()
                sat = odom[0]
                outline = time+"   "+sat+"   "+line[3:]
                wr.write(outline+"\n")
            f.close()
    wr.close()

def mjd2time(mjd):
    t0 = datetime.datetime(1858,11,17,0,0,0,0)
    return t0+datetime.timedelta(days=mjd)

## main
if len(sys.argv) != 3:
    print('#usage  : plotres.py res_filename prn')
    print('#example: plotres.py res_2021149_rov1 G01')
    sys.exit(0)

resflname = sys.argv[1]
cprn = sys.argv[2]

if (os.path.isfile('res_tmp') == False):
    if (os.path.isfile(resflname) == True):
        read_resfil(resflname)
    else:
        print('No such file:'+resflname)
        sys.exit(1)

prns = ['G01','G02','G03','G04','G05','G06','G07','G08','G09','G10','G11','G12','G13','G14','G15','G16',
      'G17','G18','G19','G20','G21','G22','G23','G24','G25','G26','G27','G28','G29','G30','G31','G32',
      'R01','R02','R03','R04','R05','R06','R07','R08','R09','R10','R11','R12','R13','R14','R15','R16',
      'R17','R18','R19','R20','R21','R22','R23','R24','E01','E02','E03','E04','E05','E06','E07','E08',
      'E09','E10','E11','E12','E13','E14','E15','E16','E17','E18','E19','E20','E21','E22','E23','E24',
      'E25','E26','E27','E28','E29','E30','E31','E32','E33','E34','E35','E36','C01','C02','C03','C04',
      'C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20',
      'C21','C22','C23','C24','C25','C26','C27','C28','C29','C30','C31','C32','C33','C34','C35','C36',
      'C37','C38','C39','C40','C41','C42','C43','C44','C45','C46','C47','C48','C56','C57','C58','C59',
      'C60','C61','J01','J02','J03','J07']
try:
    prn = prns.index(cprn)
except ValueError:
    print('Incorrect PRN format: ' + cprn)
    sys.exit(1)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    data = np.loadtxt('res_tmp', dtype=str, comments='#')

datatmp = data[data[:,2]==cprn,:]
if (len(datatmp) == 0):
    print('No satellite residuals: '+cprn)
    sys.exit(1)

print('Plotting ' + cprn)

n = len(datatmp[:,0])
dates = [0]*n
for i in range(n):
    dates[i] = mjd2time(int(datatmp[i,0]))+\
    datetime.timedelta(seconds=float(datatmp[i,1]))

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12,10))
fig.subplots_adjust(hspace=0.05)
fig.align_ylabels()

n = len(datatmp[:,0])
res = np.zeros((n,3))
for i in range(n):
    res[i,0] = float(datatmp[i,3])
    res[i,1] = float(datatmp[i,4])
    res[i,2] = float(datatmp[i,8])

ax2 = [0]*2
for i in range(2):
    axs[i].scatter(dates[:], res[:,i], marker='.',s=20, color='k')
    ax2[i] = axs[i].twinx()
    ax2[i].scatter(dates[:], res[:,2], marker='.', s=20, color='r')

axs[0].text(0.98, 0.94, cprn, horizontalalignment='right', verticalalignment='top', 
fontsize=16, transform=axs[0].transAxes)

fwidth = 1.5
for ax in axs:
    ax.spines['bottom'].set_linewidth(fwidth)
    ax.spines['left'].set_linewidth(fwidth)
    ax.spines['top'].set_linewidth(fwidth)
    ax.spines['right'].set_linewidth(fwidth)

nall = len(data[:,0])
beg = mjd2time(int(data[0,0]))+datetime.timedelta(seconds=float(data[0,1]))
end = mjd2time(int(data[nall-1,0]))+datetime.timedelta(seconds=float(data[nall-1,1]))
dt = (end-beg).total_seconds()
mjd = np.unique(data[:,0])

if dt < 86400:
    mint_intv = int(15*(math.ceil(dt/3600/3)))
    axs[1].set_xlim(beg, end)
    if mint_intv >=60:
       hour_intv=int(math.ceil((mint_intv)/60))
       print(hour_intv)
       intvs = mdates.HourLocator(interval=hour_intv)
       for ax in axs:
           ax.xaxis.set_major_locator(intvs)
else:
    nintv = int(12/len(mjd))
    hour_intv = int(math.ceil(24/nintv))
    if hour_intv%2 != 0:
        hour_intv = hour_intv+1
    intvs = mdates.HourLocator(interval=hour_intv)
    axs[1].set_xlim(beg,end)
    for ax in axs:
        ax.xaxis.set_major_locator(intvs)

for ax in axs:
    ax.tick_params(axis='both', which='both', labelsize=16, pad=4)
    ax.tick_params(axis='both', which='major', length=6)
    ax.tick_params(axis='both', which='minor', length=4)
    ax.tick_params(axis='both', which='both', direction='in', right=False, top=True)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)

## multiple days
if len(mjd) > 1:
    days=mdates.DayLocator(interval=1)
    newax=axs[1].twiny()
    newax.set_xlim(beg,end)
    newax.set_frame_on(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.tick_params(axis='both', which='major', direction='in', length=8)
    newax.tick_params(axis='both', which='minor', direction='in', length=0)
    newax.tick_params(axis='both', which='both', labelsize=16, pad=20)
    newax.xaxis.set_major_locator(days)
    newax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

min_elev = math.floor(np.nanmin(res[:,2]))
max_elev = math.ceil(np.nanmax(res[:,2]))
for ax in ax2:
    ax.tick_params(axis='y', labelsize=16, pad=4)
    ax.tick_params(axis='y', which='major', length=6)
    ax.tick_params(axis='y', which='minor', length=4)
    ax.tick_params(axis='y', which='both', direction='in', left=False)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.label.set_size(16)
    ax.set_ylim([min_elev,max_elev])

axs[0].set_ylim([-0.1, 0.1])
axs[0].yaxis.set_major_locator(MultipleLocator(0.05))
axs[0].yaxis.set_minor_locator(MultipleLocator(0.01))
axs[1].set_ylim([-5.20, 5.20])
axs[1].yaxis.set_major_locator(MultipleLocator(2))
axs[1].yaxis.set_minor_locator(MultipleLocator(0.4))

xlabel = "Time"
ylabels = ('Carrier-phase residuals (m)', 'Pseudorange residuals (m)')
for i in range(2):
    axs[i].set_ylabel(ylabels[i], labelpad=3)
    ax2[i].set_ylabel('Elevation angles (deg)', labelpad=3)
if len(mjd) > 1:
    axs[1].set_xlabel(xlabel, labelpad=23)
else:
    axs[1].set_xlabel(xlabel, labelpad=3)

pngflname=cprn+'.png'
plt.savefig(pngflname, pad_inches=None, bbox_inches='tight')
plt.show()
os.system('rm res_tmp')
sys.exit(0)
