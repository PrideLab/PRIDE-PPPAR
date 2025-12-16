#!/usr/bin/env python3

###############################################################################
##                                                                           ##
##  PURPOSE: Convert KIN file to CSV to allow loading it in QGIS             ##
##                                                                           ##
##  AUTHOR : LIF - CEDEX (derived from plottrack.py)                         ##
##                                                                           ##
##  VERSION: ver 1.0                                                         ##
##                                                                           ##
##  DATE   : Oct-10, 2025                                                    ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##    Copyright (C) 2025 by CEDEX                                            ##
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
import pandas as pd

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
  print('#usage  : plottrack.py kin_filename csv_filename')
  print('#example: plottrack.py kin_2021149_rov1 rov1_2021149')
  sys.exit(0)

kinflname = sys.argv[1]
outputflname = sys.argv[2]

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

df = pd.DataFrame(datatmp)
df  = df.rename(columns={0: "Station", 1: "* Mjd", 2: "Sod", 3: "X", 4: "Y", 5: "Z", 6: "Latitude", 7: "Longitude", 8: "Height", 9 : "Nsat/GREC2C3J_0", 10 : "Nsat/GREC2C3J_1", 11 : "Nsat/GREC2C3J_2", 12 : "Nsat/GREC2C3J_3", 13 : "Nsat/GREC2C3J_4", 14 : "Nsat/GREC2C3J_5", 15 : "Nsat/GREC2C3J_6",  16 : "PDOP"})
df.to_csv(outputflname + ".csv", index_label="ID")
os.remove('kin_tmp')
sys.exit(0)
