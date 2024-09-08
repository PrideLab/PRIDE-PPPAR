!
!! const.h  
!!
!!    Copyright (C) 2023 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program. If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
!! 
!!    Constant Parameters for PRIDE PPP-AR v3.0
!

!! ----------------------------------------- !!
!!                                           !!
!! I. Model Constant                         !!
!!                                           !!
!! ----------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! Ia. Mathematical and Physical Constant
  !!
  !! --------------------------------------------------------------------------------------------------- !!

real*8, parameter :: PI     = 3.14159265400d0
real*8, parameter :: VLIGHT = 299792458.000d0
real*8, parameter :: GM     = 398600.441800d0
real*8, parameter :: GMS    = 1.32712442076d11  !! km3/s-2
real*8, parameter :: GMM    = 4.90280070054d3   !! km3/s-2
real*8, parameter :: ERAD   = 6378.13660000d0   !! km
real*8, parameter :: OMEGA  = 7.29211500000d-5  !! rad/s

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! Ib. Time System Constant                   
  !!
  !! --------------------------------------------------------------------------------------------------- !!

real*8, parameter :: GPSTAI = 19.0000000000d0   !! s
real*8, parameter :: TAITDT = 32.1840000000d0   !! s
real*8, parameter :: GPSTDT = GPSTAI + TAITDT   !! s
real*8, parameter :: MJD2JD = 2400000.50000d0   !! d

!! ----------------------------------------- !!
!!                                           !!
!! II. General Information                   !!
!!                                           !!
!! ----------------------------------------- !!

integer*4, parameter :: MAXSYS =  5
integer*4, parameter :: MAXFRQ = 10
integer*4, parameter :: MAXATR =  9
integer*4, parameter :: MAXTYP = MAXATR * 4

character*5, parameter :: GNSS_PRIO = 'GRECJ'
character*3, parameter :: GNSS_NAME_LEN3(MAXSYS) = &
    ["GPS",     "GLO",     "GAL",     "BDS",     "QZS"    ]
character*7, parameter :: GNSS_NAME_LEN7(MAXSYS) = &
    ["GPS    ", "GLONASS", "Galileo", "BeiDou ", "QZSS   "]

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! IIa. GNSS Signal Frequency (unit: Hz)
  !! 
  !!  Sort by the frequency number in the RINEX convention
  !!
  !! ---- L1 -------- L2 -------------------- L5 -------- L6 -------- L7 --------- L8 ------------------ !!

real*8, parameter :: FREQ_G(MAXFRQ) = &
    [1.575420d9, 1.227600d9, 0.d0, 0.d0, 1.176450d9,       0.d0,       0.d0,       0.d0, 0.d0, 0.d0]
real*8, parameter :: FREQ_R(MAXFRQ) = &
    [1.602000d9, 1.246000d9, 0.d0, 0.d0,       0.d0,       0.d0,       0.d0,       0.d0, 0.d0, 0.d0]
real*8, parameter :: FREQ_E(MAXFRQ) = &
    [1.575420d9,       0.d0, 0.d0, 0.d0, 1.176450d9, 1.278750d9, 1.207140d9, 1.191795d9, 0.d0, 0.d0]
real*8, parameter :: FREQ_C(MAXFRQ) = &
    [1.575420d9, 1.561098d9, 0.d0, 0.d0, 1.176450d9, 1.268520d9, 1.207140d9, 1.191795d9, 0.d0, 0.d0]
real*8, parameter :: FREQ_J(MAXFRQ) = &
    [1.575420d9, 1.227600d9, 0.d0, 0.d0, 1.176450d9, 1.278750d9,       0.d0,       0.d0, 0.d0, 0.d0]
real*8, parameter :: FREQ_SYS(MAXFRQ, MAXSYS) = &
    reshape([FREQ_G, FREQ_R, FREQ_E, FREQ_C, FREQ_J], [MAXFRQ, MAXSYS])

character*3, parameter :: FREQ_NAME_G(MAXFRQ) = &
    ["L1 ",      "L2 ",  "   ",  "   ",  "L5 ",      "   ",      "   ",      "   ",    "   ",  "   "]
character*3, parameter :: FREQ_NAME_R(MAXFRQ) = &
    ["G1 ",      "G2 ",  "   ",  "   ",  "   ",      "   ",      "   ",      "   ",    "   ",  "   "]
character*3, parameter :: FREQ_NAME_E(MAXFRQ) = &
    ["E1 ",      "   ",  "   ",  "   ",  "E5a",      "E6 ",      "E5b",      "E5 ",    "   ",  "   "]
character*3, parameter :: FREQ_NAME_C(MAXFRQ) = &
    ["B1C",      "B1I",  "   ",  "   ",  "B2a",      "B3I",      "B2b",      "B2 ",    "   ",  "   "]
character*3, parameter :: FREQ_NAME_J(MAXFRQ) = &
    ["L1 ",      "L2 ",  "   ",  "   ",  "L5 ",      "L6 ",      "   ",      "   ",    "   ",  "   "]
character*3, parameter :: FREQ_NAME_SYS(MAXFRQ, MAXSYS) = &
    reshape([FREQ_NAME_G, FREQ_NAME_R, FREQ_NAME_E, FREQ_NAME_C, FREQ_NAME_J], [MAXFRQ, MAXSYS])

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! IIb. GNSS Observation Attribute (tracking mode or channel) 
  !!
  !! --------------------------------------------------------------------------------------------------- !!

character*9, parameter :: OBS_PRIO_G = 'XSLCPWQI'
character*9, parameter :: OBS_PRIO_R = 'CP'
character*9, parameter :: OBS_PRIO_E = 'CQX'
character*9, parameter :: OBS_PRIO_C = 'PXZCDI'
character*9, parameter :: OBS_PRIO_J = 'XZSLCQI'
character*9, parameter :: OBS_PRIO_SYS(MAXSYS) = [OBS_PRIO_G, OBS_PRIO_R, OBS_PRIO_E, OBS_PRIO_C, OBS_PRIO_J]

!! ----------------------------------------- !!
!!                                           !!
!! III. Processing Setting                   !!
!!                                           !!
!! ----------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! IIIa. Maximum Satellite Number 
  !!
  !! --------------------------------------------------------------------------------------------------- !!

integer*4, parameter :: MAXSAT_G = 40
integer*4, parameter :: MAXSAT_R = 30
integer*4, parameter :: MAXSAT_E = 40
integer*4, parameter :: MAXSAT_C = 70
integer*4, parameter :: MAXSAT_J = 10
integer*4, parameter :: MAXSAT_SYS(MAXSYS) = [MAXSAT_G, MAXSAT_R, MAXSAT_E, MAXSAT_C, MAXSAT_J]
integer*4, parameter :: MAXSAT = MAXSAT_G + MAXSAT_R + MAXSAT_E + MAXSAT_C + MAXSAT_J

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! IIIb. Allocation Limit and Maximum Parameter Number
  !!
  !! --------------------------------------------------------------------------------------------------- !!

integer*4, parameter :: MAXDAY = 108
integer*4, parameter :: MAXEPH = MAXDAY * (MAXSAT*24 + MAXSAT_E*120)
integer*4, parameter :: MAXPAR_STA = 5 * MAXDAY
integer*4, parameter :: MAXOW_ST = 5 * MAXDAY * MAXSAT
integer*4, parameter :: MAXPAR = MAXPAR_STA + MAXOW_ST
integer*4, parameter :: MAXSD_SIT = 40 * MAXDAY * MAXSAT

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! IIIc. Geogrephic Grid Limit
  !!
  !! --------------------------------------------------------------------------------------------------- !!

integer*4, parameter :: MAXLAT = 71
integer*4, parameter :: MAXLON = 73

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! IIId. Time Window Limit
  !!
  !! --------------------------------------------------------------------------------------------------- !!

real*8, parameter :: MAXWND = 0.01d0

  !! --------------------------------------------------------------------------------------------------- !!
  !!
  !! IIIe. Oceanload Setting 
  !!
  !! --------------------------------------------------------------------------------------------------- !!

integer*4, parameter :: MAXDEG = 100