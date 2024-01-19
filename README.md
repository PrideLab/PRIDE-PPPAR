![pridelab.icon](https://github.com/PrideLab/PRIDE-PPPAR/blob/master/doc/PRIDE.png)

## PRIDE-PPPAR ver. 3.0 (last updated on 2024-01-19)

PRIDE PPP-AR originates in Dr. Maorong Ge’s efforts on PPP-AR and later developed and improved by Dr. Jianghui Geng's team. It is an open-source software package which is based on many GNSS professionals’ collective work in GNSS Research Center, Wuhan University. We would like to thank them all for their brilliant contributions to this software. 

We make this package open source with the goal of benefiting those professionals in their early career, and also advocate the geodetic and geophysical applications of GNSS PPP-AR. Especially, we hope that this package can contribute to high-precision applications in geosciences such as crustal motion and troposphere sounding studies. The entire open-source project is funded by the National Natural Science Foundation of China (No. 42025401) and is under the auspices of IAG SC 4.2 "Positioning and Applications".

The GNSS products are accessible at <ftp://igs.gnsswhu.cn/pub/whu/phasebias>.
Latest updates for Support, Training courses and FAQ can be found at <http://pride.whu.edu.cn>.
The copyright of this package is protected by GNU General Public License (version 3). 

Relevant publications are

* Geng J, Zhang Q, Li G, Liu J, Liu D (2022). Observable-specific phase biases of Wuhan multi-GNSS experiment analysis center's rapid satellite products. *Satellite Navigation*, 3(1):1-15. doi:[10.1186/s43020-022-00084-0](https://doi.org/10.1186/s43020-022-00084-0)
* Geng J, Wen Q, Zhang Q, Li G, Zhang K (2022). GNSS observable-specific phase biases for all-frequency PPP ambiguity resolution. *Journal of Geodesy*, 96(11):1-18. doi:[10.1007/s00190-022-01602-3](https://doi.org/10.1007/s00190-022-01602-3)
* Geng J, Mao S (2021). Massive GNSS network analysis without baselines: Undifferenced ambiguity resolution. *Journal of Geophysical Research: Solid Earth*. 126(10), e2020JB021558. doi:[10.1029/2020JB021558](https://doi.org/10.1029/2020JB021558)
* Geng J, Yang S, Guo J (2021). Assessing IGS GPS/Galileo/BDS-2/BDS-3 phase bias products with PRIDE PPP-AR. *Satellite Navigation*, 2(1):1-15. doi:[10.1186/s43020-021-00049-9](https://doi.org/10.1186/s43020-021-00049-9)
* Geng J, Chen X, Pan Y, Zhao Q (2019). A modified phase clock/bias model to improve PPP ambiguity resolution at Wuhan University. *Journal of Geodesy*, 93(10):2053-2067. doi:[10.1007/s00190-019-01301-6](https://doi.org/10.1007/s00190-019-01301-6)
* Geng J, Chen X, Pan Y, Mao S, Li C, Zhou J, Zhang K (2019). PRIDE PPP‑AR: an open‑source software for GPS PPP ambiguity resolution. *GPS Solutions*, 23(91):1-10. doi:[10.1007/s10291-019-0888-1](https://doi.org/10.1007/s10291-019-0888-1)

PRIDE PPP-AR ver. 3.0 is available for:

1) Multi-GNSS data processing with GPS, GLONASS, Galileo, BDS-2/3 and QZSS
2) All-frequency PPP-AR on any dual-frequency ionosphere-free combinations of GPS, Galileo, and BDS-2/3
3) High-rate GNSS data processing up to 50 Hz
4) High-dynamic mobile platforms applicable for aerial photogrammetry, ship-borne gravimetry, etc.
5) Kinematic orbiting for LEO satellites
6) Multi-day processing up to 108 days (at a sampling rate of 30 s) without day-boundary discontinuities
7) IPPP clock estimation for time and frequency transfer
8) Troposphere modeling with Vienna Mapping Function 1/3 (VMF1/VMF3)
9) Second-order ionospheric delay correction with GIM products
10) Receiver clock jump detection and mitigation
11) Adoption of the latest IGS conventions: Bias-SINEX, IGS20 reference frame, ORBEX, RINEX 4, etc.
12) User-friendly operation and visualization for early-career researchers with lite-version GUI
13) Ambiguity-float PPP with backward compatibility from 1994 when selective availablity (SA) was on
14) Timely data processing with rapid (RAP) products and real-time archived (RTS) products

The improvements made in PRIDE PPP-AR version 3.0 include:

* Enable all-frequency PPP-AR on any dual-frequency ionosphere-free combinations
* Employ the latest WUM0MGXRAP products to resolve ambiguities on new GNSS signals (L5/E6/E5b/E5B/B1C/B2a/B2)
* Support kinematic orbiting for LEO satellites
* Improve capability and long-term consistency of multi-day processing
* Provide more command-line options and models for parameter estimation
* Increase compatibility with the latest IGS data and product extensions
* Refine the data format of result files and the output information of program runs

`Notes: The all-frequency satellite orbit/clock/bias/ERP/quaternion products from Wuhan University are required by PRIDE PPP-AR`

## Version History

See our [Change Log](https://github.com/PrideLab/PRIDE-PPPAR/blob/master/CHANGELOG.md) for detailed update history before version 3.0.

### 2023-12-13 (v3.0)

* `pdp3`: add IGN FTP site <ftp://igs.ign.fr/pub/igs/products/mgex/> to GNSS product download path
* `pdp3`: refine the procedures for ambiguity resolution when missing code/phase biases for some satellites

### 2023-11-24 (v3.0)

Release of **PRIDE PPP-AR v3.0**

* Enable **all-frequency PPP-AR** on any dual-frequency ionosphere-free combinations
  * `pdp3`: add `-frq` option to choose the signal frequencies for each GNSS
  * `pdp3`: add verification for the availability of observations and code/phase OSBs
* Employ the latest **WUM0MGXRAP** products to resolve ambiguities on new GNSS signals
  * Rapid all-frequency phase clock/bias products from Wuhan Univeristy
  * Support ambiguity resolution on GPS L1/L2/L5, Galileo E1/E5a/E6/E5b/E5, BDS B1C/B1I/B2a/B3I/B2
  * Append satellite orbit, clock, and quaternion records at the second midnight
  * Align phase clock/bias across consecutive days to mitigate day-bounday discontinuities
* Support precise orbit determination for **LEO satellites**
  * `pdp3`: add `L` positioning mode for kinematic orbiting of LEO satellites
  * `table`: add processing scripts for GRACE and GRACE-FO, other LEO satellites can be integrated manually
  * `table`: add setting `LEO quaternions` in `config_template` to specify attitude products for LEO satellites
* Improve capability and long-term consistency of **multi-day processing**
  * Enhance capability of multi-day processing to 108 days
  * Enable establishing integer-ambiguity constraints across midnights with aligned GNSS products to **mitigate day-boundary discontinuities** in the result
  * `table`: add setting `Truncate at midnight` in `config_template` to set whether to truncate ambiguities at the day boundaries
* Provide more command-line options and models for parameter estimation
  * `pdp3`: add `-r` option to choose `WNO`/`STO` model for receiver clock estimation
  * `pdp3`: add `-h` option to choose `NON`/`PWC`/`STO` model for HTG estimation
* Increase compatibility with the latest IGS data and product extensions
  * Support observation data and broadcast ephemerides in **RINEX 4**
  * `lib`: fix bug for not reading interspersing AR and AS clock records
  * `lib`: fix bug for not reading time-varying OSB estimates and their gradients
  * `lib`: fix bug for not reading non-16-digit quaternions
* Refine the data format of result files and the output information of program runs
  * Add header information (antenna, PCO, frequency selection, positioning mode, etc.) in result files
  * Increase the number of digits in `kin` file and `res` file
  * `pdp3`: add `-v` option to output detailed information from `arsig` and `lsq`
  * `xyz2enu`: output RMS instead of STD when reference coordinates are specified

Special notes for compatibility with older versions of PRIDE PPP-AR:

1) `pdp3_Mac.sh` has been integrated into `pdp3.sh`
2) Older versions of PRIDE PPP-AR cannot recognize the new namings of WUM0MGXRAP and WUM0MGXRTS products
3) Older versions of PRIDE PPP-AR cannot read the new settings in configuration file
4) Older versions of PRIDE PPP-AR cannot read the new format in result files
 
## Getting in Touch

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FPrideLab%2FPRIDE-PPPAR&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=PAGE+VIEWS&edge_flat=false)](https://hits.seeyoufarm.com)

* You can contact us for **bug reports** and **comments** by sending an email or leaving a message on our website:
  * Email: <pride@whu.edu.cn>
  * Website: <http://pride.whu.edu.cn>
* For Chinese users, we provide Tencent **QQ Group** service.
  * QQ group: 971523302

## License

***Copyright (C) 2023 by Wuhan University, All rights reserved.***
