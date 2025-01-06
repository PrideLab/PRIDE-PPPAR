#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Prepare LEO RINEX observation files and products                ##
##                                                                           ##
##  AUTHOR : the PRIDE Group pride@whu.edu.cn                                ##
##                                                                           ##
##  VERSION: ver 3.1.0                                                        ##
##                                                                           ##
##  DATE   : Jan-04, 2025                                                    ##
##                                                                           ##
##              @ GNSS RESEARCH CENTER, WUHAN UNIVERSITY, 2025               ##
##                                                                           ##
##    Copyright (C) 2025 by Wuhan University                                 ##
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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


## main
if len(sys.argv) != 2 and len(sys.argv) != 3:
    print('#usage  : plotmhm.py mhm_filename png_filename')
    print('#example: plotmhm.py mhm_rov1 png_mhm_rov1')
    print('#if no png_filename, default mhm_filename_P or mhm_filename_L')
    sys.exit(0)

mhmflname = sys.argv[1]
if len(sys.argv) == 3:
    pngflname = sys.argv[2]
elif len(sys.argv) == 2:
    pngflname = sys.argv[1]

if os.path.isfile(mhmflname) == True:
    with open(mhmflname, 'r') as f:
        ## read header
        site = "NONE"
        daystart = None
        dayend = None
        while (1):
            line = f.readline().strip("\n")
            if "STATION" in line:
                site = line[0:4]
            elif "DAY LIST" in line:
                doys = line.split()
                if daystart is None:
                    daystart = doys[0]
                dayend = doys[-3]
            elif " MODELING DURATION (DAY)" in line:
                daynums = line.split()
                daynum = daynums[0]
            elif "END OF HEADER" in line:
                break
        ## prepare data
        sysnum = 0
        sysname = []
        dataname = []
        if line == "":
            f.close()
        else:
            while (1):
                line = f.readline().strip("\n")
                if line == "":
                    break
                elif line[0] == "*":
                    continue
                elif "START OF" in line:
                    commonnum = 0
                    line = line.split()
                    for i in range(2, len(line), 2):
                        if line[i] != "QZS":
                            sysname.append(line[i])
                            dataname.append('data_' + str(sysnum))
                            commonnum += 1
                            sysnum += 1
                    datatmp = []
                    f.readline()
                    while (1):
                        line = f.readline().strip("\n")
                        if "END OF" in line:
                            break
                        else:
                            line = line.split()
                            AZ = float(line[0]) / 180 * math.pi
                            EL = 90 - float(line[1])
                            PCOR = float(line[2])
                            LCOR = float(line[6])
                            dataline = [AZ, EL, PCOR, LCOR]
                            datatmp.append(dataline)
                    for i in range(sysnum-commonnum, sysnum):
                        locals()[dataname[i]] = datatmp
    f.close()
else:
    print('No such file:' + mhmflname)
    sys.exit(1)

del datatmp

# title
axtitle = []
for i in range(sysnum):
    if sysname[i] == 'GPS':
        axtitle.append('GPS')
    elif sysname[i] == 'GAL':
        axtitle.append('Galileo')
    elif sysname[i] == 'BDS':
        axtitle.append('BDS')
    elif sysname[i] == 'GLO':
        axtitle.append('Glonass')
    else:
        axtitle.append(sysname[i])

# parameters for colorbar
if sysnum == 1:
    figwidth = 6
    figheight = 5
    asp = 40
    shr = 0.9
elif sysnum == 2:
    figwidth = 6
    figheight = 5
    asp = 34
    shr = 0.75
elif sysnum == 3:
    figwidth = 9
    figheight = 5
    asp = 42
    shr = 0.75
elif sysnum == 4:
    figwidth = 12
    figheight = 5
    asp = 44
    shr = 0.75
else:
    figwidth = 3 * sysnum
    figheight = 5
    asp = 42
    shr = 0.75

# plotmhm Pseudorange
fig1, axs = plt.subplots(1, sysnum, subplot_kw=dict(projection='polar'), figsize=(figwidth, figheight))
fig1.subplots_adjust(wspace=0.25, hspace=0)
if sysnum == 1:
    datatmp = np.array(locals()[dataname[0]])
    scatter = axs.scatter(datatmp[:, 0], datatmp[:, 1], c=datatmp[:, 2], s=6, cmap='jet', vmin=-1.5, vmax=1.5, alpha=1)
    del datatmp
    axs.set_title(axtitle[0], fontsize=18, pad=34)
    axs.set_theta_zero_location('N')
    axs.set_theta_direction(-1)
    axs.set_thetagrids(np.arange(0, 360, 60), fontsize=10)
    axs.set_rgrids(np.arange(0, 90, 30))
    axs.set_rlabel_position(0)
    axs.set_rlim(0, 90)
    axs.set_yticklabels(['', '60\xb0', '30\xb0'], fontsize=10)
    axs.grid(True)
elif sysnum > 1:
    for i in range(sysnum):
        datatmp = np.array(locals()[dataname[i]])
        scatter = axs[i].scatter(datatmp[:, 0], datatmp[:, 1], c=datatmp[:, 2], s=6, cmap='jet', vmin=-1.5, vmax=1.5, alpha=1)
        del datatmp
    for i in range(sysnum):
        axs[i].set_title(axtitle[i], fontsize=18, pad=34)
        axs[i].set_theta_zero_location('N')  # 将极轴的起始位置设置为正北方向
        axs[i].set_theta_direction(-1)
        axs[i].set_thetagrids(np.arange(0, 360, 60), fontsize=10)
        axs[i].set_rgrids(np.arange(0, 90, 30))
        axs[i].set_rlabel_position(0)
        axs[i].set_rlim(0, 90)
        axs[i].set_yticklabels(['', '60\xb0', '30\xb0'], fontsize=10)
        axs[i].grid(True)

# colorbar
cbar = fig1.colorbar(scatter, ax=axs, orientation='horizontal', aspect=asp, shrink=shr, pad=0.12)
cbar.ax.xaxis.set_tick_params(labelsize=10)
norm = Normalize(vmin=-1.5, vmax=1.5)
cbar.set_ticks([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])

site = site.upper()
if int(daynum) == 1:
    titlepse = 'Pseudorange Multipath Correction (meter)\n' + site + ' ' + daystart
elif int(daynum) > 1:
    titlepse = 'Pseudorange Multipath Correction (meter)\n' + site + ' ' + daystart + '-' + dayend

# title
plt.suptitle(titlepse, fontsize=20, y=0.12)
plt.savefig(pngflname + '_P', bbox_inches='tight', dpi=250)


# plotmhm Carrier Phase
fig2, axs = plt.subplots(1, sysnum, subplot_kw=dict(projection='polar'), figsize=(figwidth, figheight))
fig2.subplots_adjust(wspace=0.25, hspace=0)
if sysnum == 1:
    datatmp = np.array(locals()[dataname[0]])
    scatter = axs.scatter(datatmp[:, 0], datatmp[:, 1], c=datatmp[:, 3], s=6, cmap='jet', vmin=-0.03, vmax=0.03, alpha=1)
    del datatmp
    axs.set_title(axtitle[0], fontsize=18, pad=34)
    axs.set_theta_zero_location('N')
    axs.set_theta_direction(-1)
    axs.set_thetagrids(np.arange(0, 360, 60), fontsize=10)
    axs.set_rgrids(np.arange(0, 90, 30))
    axs.set_rlabel_position(0)
    axs.set_rlim(0, 90)
    axs.set_yticklabels(['', '60\xb0', '30\xb0'], fontsize=10)
    axs.grid(True)
elif sysnum > 1:
    for i in range(sysnum):
        datatmp = np.array(locals()[dataname[i]])
        scatter = axs[i].scatter(datatmp[:, 0], datatmp[:, 1], c=datatmp[:, 3], s=6, cmap='jet', vmin=-0.03, vmax=0.03, alpha=1)
        del datatmp
        axs[i].set_title(axtitle[i], fontsize=18, pad=34)
        axs[i].set_theta_zero_location('N')  # 将极轴的起始位置设置为正北方向
        axs[i].set_theta_direction(-1)
        axs[i].set_thetagrids(np.arange(0, 360, 60), fontsize=10)
        axs[i].set_rgrids(np.arange(0, 90, 30))
        axs[i].set_rlabel_position(0)
        axs[i].set_rlim(0, 90)
        axs[i].set_yticklabels(['', '60\xb0', '30\xb0'], fontsize=10)
        axs[i].grid(True)

# colorbar
cbar = fig2.colorbar(scatter, ax=axs, orientation='horizontal', aspect=asp, shrink=shr, pad=0.12)
cbar.ax.xaxis.set_tick_params(labelsize=10)
norm = Normalize(vmin=-0.1, vmax=0.1)
cbar.set_ticks([-0.03, -0.015, 0, 0.015, 0.03])

site = site.upper()
if int(daynum) == 1:
    titlecar = 'Carrier Phase Multipath Correction (meter)\n' + site + ' ' + daystart
elif int(daynum) > 1:
    titlecar = 'Carrier Phase Multipath Correction (meter)\n' + site + ' ' + daystart + '-' + dayend

# title
plt.suptitle(titlecar, fontsize=20, y=0.12)
plt.savefig(pngflname + '_L', bbox_inches='tight', dpi=250)

sys.exit(0)

