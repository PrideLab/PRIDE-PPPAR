#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Plot kinematic position of LEO satellite                        ##
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
import math
from datetime import datetime, timedelta
import platform

try:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.ticker import FormatStrFormatter, MultipleLocator
except ModuleNotFoundError:
    print("Error: %s need 'numpy' and 'matplotlib' modules of python for plotting!" % sys.argv[0])
    sys.exit(1)

global PATHSEP
if platform.system() == "Windows":
    PATHSEP = "\\"
else:
    PATHSEP = '/'


# Time function
def mjd2ydoy(mjd):
    time = mjd2time(mjd)
    time0 = datetime(time.year, 1, 1, 0, 0, 0)
    return time.year, (time - time0).days + 1


def time2mjd(t):  # t=[y, m, d, h, min, sec]
    time0 = datetime(1858, 11, 17, 0, 0, 0, 0)
    time = datetime(t[0], t[1], t[2], 0, 0, 0)
    mjd = (time-time0).days
    sod = t[3]*3600+t[4]*60+t[5]
    return mjd, sod


def mjd2time(mjd):
    time0 = datetime(1858, 11, 17, 0, 0, 0, 0)
    time1 = time0 + timedelta(days=mjd)
    return time1


def mjdsod2time(mjd, sod):
    time0 = datetime(1858, 11, 17, 0, 0, 0, 0)
    time1 = time0 + timedelta(days=mjd)
    time1 = time1 + timedelta(seconds=int(sod))
    fsod = sod - int(sod)
    millsecond = int(fsod*1e3)
    time1 = time1 + timedelta(milliseconds=millsecond)
    return time1


def gfot_mjd(grasec):  # GRACE/GRACE-FO Time stamp
    time0 = datetime(2000, 1, 1, 12, 00, 00)
    time1 = time0 + timedelta(seconds=grasec)
    fsec = grasec-int(grasec)
    t = [time1.year, time1.month, time1.day, time1.hour, time1.minute, time1.second+fsec]
    mjd, sod = time2mjd(t)
    return mjd, sod


# def gpst2utc():  # UTC+leap.sec=TAI GPST+19=TAI
# Vector operation
def rmse(a):
    rms = 0
    numzero = 0
    for i in a:
        if np.isnan(i):
            numzero += 1
            continue
        rms = rms + pow(i, 2)
        if i == 0:
            numzero += 1
    if numzero == len(a):
        rms = np.nan
    else:
        rms = np.sqrt(rms/(len(a)-numzero))
    return rms


def cross(a, b):  # 3d vector cross product
    c = np.zeros(3)
    c[0] = a[1]*b[2]-a[2]*b[1]
    c[1] = a[2]*b[0]-a[0]*b[2]
    c[2] = a[0]*b[1]-a[1]*b[0]
    return c


def unit_vector(a, n):
    length = 0
    for i in range(n):
        length = length + a[i]*a[i]
    length = np.sqrt(length)
    if length == 0:
        sys.exit()
    return a/length


def flitter(a, ratio):  # a(n, 3)
    b = a
    num = 0
    for i in range(len(a)):
        for j in range(3):
            if abs(a[i, j]) > ratio[j]:
                b[i, :] = np.nan
        if np.isnan(b[i, 0]) and np.isnan(b[i, 1]) and np.isnan(b[i, 2]):
            num += 1
    return b, num/len(b)


# PSO and kin
class ORBIT:
    def __init__(self):
        self.mjd = 0
        self.sod = 0.0
        self.xyz = np.zeros(3)
        self.vel = np.zeros(3)
        self.satnum = 0
        self.pdop = 0
        self.mwin = 0.5  # time match window

    def time_match(self, mjd, sod):
        lmatch = False
        if self.mjd == mjd and abs(self.sod-sod) < self.mwin:
            lmatch = True
        return lmatch

    def uv_rtn(self):
        e_a = unit_vector(self.vel, 3)
        e_c = cross(self.xyz, self.vel)
        e_c = unit_vector(e_c, 3)
        e_r = cross(e_a, e_c)
        return e_a, e_c, e_r


# Read files
def read_kin(kinfl):
    try:
        allmjd = []
        with open(kinfl, 'r') as rdkin:
            # Header
            while 1:
                line = rdkin.readline().strip('\n')
                if 'STATION' in line:
                    strs = line.split()
                    leocode = strs[0]
                if 'END OF HEADER' in line:
                    rdkin.readline()
                    break
            # Data set
            datatmp = rdkin.readlines()
        num = len(datatmp)
        kin_data = [ORBIT() for i in range(num)]
        for i in range(num):
            strs = datatmp[i].split()
            if '*' in datatmp[i]:
                kin_data[i].mjd = int(strs[0])
                kin_data[i].sod = float(strs[1])
                kin_data[i].xyz = np.nan
                kin_data[i].satnum = np.nan
                kin_data[i].pdop = np.nan
                continue
            kin_data[i].mjd = int(strs[0])
            kin_data[i].sod = float(strs[1])
            kin_data[i].xyz = [float(strs[j+2]) for j in range(3)]
            kin_data[i].satnum = float(strs[8])
            kin_data[i].pdop = float(strs[15])
            if kin_data[i].mjd not in allmjd:
                allmjd.append(kin_data[i].mjd)
        return kin_data, leocode, allmjd
    except IOError:
        print("Error : open kin_ file \"" + os.path.abspath(kinfl)+"\" failed!")
        sys.exit(1)


def read_grapso(rdgra, lheader, kmjd, ksod):  # GRACE PSO: GNV
    if lheader:
        while 1:
            line = rdgra.readline().strip('\n')
            if '# End of YAML header' in line or 'END OF HEADER' in line:
                break
    else:
        pso_data = ORBIT()
        while 1:
            findex = rdgra.tell()
            line = rdgra.readline().strip('\n')
            if line == '':
                break
            strs = line.split()
            grasec = int(strs[0])
            mjd, sod = gfot_mjd(grasec)
            if mjd > kmjd or mjd == kmjd and sod-ksod > pso_data.mwin:
                rdgra.seek(findex)
                break
            pso_data.mjd = mjd
            pso_data.sod = sod
            pso_data.xyz = [float(strs[i+3]) for i in range(3)]
            pso_data.vel = [float(strs[i+9]) for i in range(3)]
            # if pso_data.time_match(kmjd, ksod):
            if kmjd == mjd and ksod == sod:
                break
        return pso_data


def read_psoh(rdpso, leocode):
    if leocode[:3] == 'gra':
        read_grapso(rdpso, True, -1, -1)


def read_psodata(rdpso, leocode, mjd, sod):
    if leocode[:3] == 'gra':
        psodata = read_grapso(rdpso, False, mjd, sod)
    return psodata


# Plot
def plot_rtn(rtn, pngfl):
    plt.rc('font', size=24)
    fig, axs = plt.subplots(4, 1, figsize=(16, 12), sharex=True)
    fig.subplots_adjust(hspace=0.3)
    colors = ('r', 'g', 'b')
    labels = ('Along', 'Corss', 'Radial')
    # ftsize = 18

    dates = []
    for i in range(len(rtn)):
        mjd = rtn[i, 0]
        sod = rtn[i, 1]
        time = mjdsod2time(mjd, sod)
        dates.append(time)
    mjds = np.unique(rtn[:, 0])
    dt = (dates[-1] - dates[0]).total_seconds()
    if len(mjds) == 1:
        mint_intv = int(15*(math.ceil(dt/3600/3)))
        intvs = mdates.MinuteLocator(interval=mint_intv)
        if mint_intv >= 60:
            hour_intv = int(math.ceil((mint_intv)/60))
            intvs = mdates.HourLocator(interval=hour_intv)
    else:
        if len(mjds) <= 12:
            nintv = int(12/len(mjds))
        else:
            nintv = 1
        hour_intv = int(math.ceil(24/nintv))
        if hour_intv % 2 != 0:
            hour_intv = hour_intv+1
        intvs = mdates.HourLocator(interval=hour_intv)

    for i in range(3):
        rms = rmse(rtn[:, i+2]*100)
        mean = np.nanmean(rtn[:, i+2]*100)
        # axs[i].scatter(rtn[:, 1], rtn[:, i+2]*100, label=("RMS %5.2fcm MEAN %5.2fcm" % (rms, mean)), color=colors[i], s=2, marker="o")
        axs[i].scatter(dates, rtn[:, i+2]*100, label=("RMS %5.2fcm MEAN %5.2fcm" % (rms, mean)), color=colors[i], s=2, marker="o")
        axs[i].set(ylabel=labels[i]+' (cm)', xlim=[dates[0], dates[-1]])  # , ylim=[np.nanmin(rtn[:, i+2]*100)-3, np.nanmax(rtn[:, i+2])*100+3])
        axs[i].legend(fontsize=16, loc=1, frameon=True, markerscale=4)  # , bbox_to_anchor=(0.6, 1.2), borderaxespad=0.)
    # left = axs[3].plot(rtn[:, 1], rtn[:, 5], color='cyan', label="Satnum")
    left = axs[3].plot(dates, rtn[:, 5], color='cyan', label="Satnum")
    axs[3].tick_params(axis='y', colors='cyan')
    axs[3].spines["left"].set_color('cyan')
    axs[3].set(ylabel='Satnum', xlim=[dates[0], dates[-1]], ylim=[0, np.nanmax(rtn[:, 5])+3])
    axr = axs[3].twinx()
    # right = axr.plot(rtn[:, 1], rtn[:, 6], color='y', label="PDOP")
    right = axr.plot(dates, rtn[:, 6], color='y', label="PDOP")
    axr.tick_params(axis='y', colors='y')
    axr.spines["right"].set_color('y')
    axr.set(ylabel='PDOP', xlim=[dates[0], dates[-1]])  # , ylim=[0, 10])
    doubley = left+right
    labels = [ax.get_label() for ax in doubley]
    axs[3].legend(doubley, labels, fontsize=16, loc=1, frameon=False, ncol=2)

    for i in range(5):
        if i == 4:
            ax = axr
        else:
            ax = axs[i]
        ax.tick_params(axis='both', which='both', labelsize=16, pad=4)
        ax.tick_params(axis='both', which='major', length=6)
        ax.tick_params(axis='both', which='minor', length=4)
        ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
        ax.xaxis.set_major_locator(intvs)
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    if len(mjds) > 1:
        days = mdates.DayLocator()
        newax = axs[3].twiny()
        newax.set_xlim(dates[0], dates[-1])
        newax.set_frame_on(False)
        newax.xaxis.set_ticks_position('bottom')
        newax.xaxis.set_label_position('bottom')
        newax.tick_params(axis='x', which='major', direction='in', length=8)
        newax.tick_params(axis='x', which='minor', direction='in', length=0)
        newax.tick_params(axis='both', which='both', labelsize=16, pad=20)
        newax.xaxis.set_major_locator(days)
        newax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    xlabel = "Time"
    if len(mjds) > 1:
        axs[3].set_xlabel(xlabel, labelpad=23)
    else:
        axs[3].set_xlabel(xlabel, labelpad=3)
    plt.savefig(pngfl, pad_inches=None, bbox_inches='tight')
    plt.show()


def open_allpso(psofl0, allmjd):
    if not os.access(psofl0, os.R_OK):
        print("Error : the input pso file \"" + "\" is unreadable!")
        sys.exit(1)

    # pso_path = os.path.abspath(psofl0)
    # pso_path = os.path.relpath(psofl0)
    pso_path = os.path.dirname(psofl0)
    leosat = psofl0.split('_')[2]

    readpsos = {}
    for mjd in allmjd:
        year, day = mjd2ydoy(mjd)
        psofl = pso_path + PATHSEP + "pso_%4d%03d_%s" % (year, day, leosat)
        try:
            readpso = open(psofl, 'r')
        except IOError:
            print("Warning : open pso file for " + str(mjd) + "(" + psofl + ") failed!")
            readpso = None
        readpsos[mjd] = readpso
    return readpsos


# Main
def main(kinfl, psofl, pngfl):
    kin_data, leocode, allmjd = read_kin(kinfl)
    rdpsos = open_allpso(psofl, allmjd)
    index_header = {}
    for mjd in allmjd:
        index_header[mjd] = False
        if rdpsos[mjd] is not None:
            index_header[mjd] = True
    rtn = np.zeros((len(kin_data), 7))
    i = 0
    for skin_data in kin_data:
        # if skin_data.mjd == 0 and skin_data.sod == 0.0:
        if np.isnan(skin_data.satnum) and np.isnan(skin_data.pdop):
            rtn[i, 0] = skin_data.mjd
            rtn[i, 1] = skin_data.sod
            rtn[i, 2:] = np.nan
            i = i + 1
            continue
        rdpso = rdpsos[skin_data.mjd]
        if rdpso is None:
            continue
        if index_header[skin_data.mjd]:  # header
            read_psoh(rdpso, leocode)
            index_header[skin_data.mjd] = False
        pso_data = read_psodata(rdpso, leocode, skin_data.mjd, skin_data.sod)
        if pso_data.mjd == 0 and pso_data.sod == 0.0:
            rtn[i, 0] = skin_data.mjd
            rtn[i, 1] = skin_data.sod
            rtn[i, 2:] = np.nan
            i = i + 1
            continue
        rtn[i, 0] = skin_data.mjd
        rtn[i, 1] = skin_data.sod
        e_a, e_c, e_r = pso_data.uv_rtn()
        xyzd = [pso_data.xyz[j]-skin_data.xyz[j] for j in range(3)]
        rtn[i, 2] = np.dot(e_a, xyzd)  # track/along
        rtn[i, 3] = np.dot(e_c, xyzd)  # normal/corss
        rtn[i, 4] = np.dot(e_r, xyzd)  # radial
        rtn[i, 5] = skin_data.satnum
        rtn[i, 6] = skin_data.pdop
        i = i + 1
    # rms_rtn = [rmse(rtn[:, i+2]) for i in range(3)]
    # print(rms_rtn[0], rms_rtn[1], rms_rtn[2])
    plot_rtn(rtn, pngfl)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print('Usage: kinfl psofl pngfl')
        sys.exit(0)
    kinfl = sys.argv[1]
    psofl = sys.argv[2]
    pngfl = sys.argv[3]
    main(kinfl, psofl, pngfl)
