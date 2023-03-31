![pridelab.icon](https://github.com/PrideLab/PRIDE-PPPAR/blob/master/doc/PRIDE.png)
## PRIDE-PPPAR ver. 2.2

PRIDE PPP-AR ver. 2.2 originates in Dr. Maorong Ge’s efforts on PPP-AR and later developed and improved by Dr. Jianghui Geng's team. It is an open-source software package which is based on many GNSS professionals’ collective work in GNSS Research Center, Wuhan University. We would like to thank them all for their brilliant contributions to this software. 

We make this package open source with the goal of benefiting those professionals in their early career, and also advocate the geodetic and geophysical applications of GNSS PPP-AR. Especially, we hope that this package can contribute to high-precision applications in geosciences such as crustal motion and troposphere sounding studies. The entire open-source project is funded by the National Natural Science Foundation of China (No. 42025401) and is under the auspices of IAG SC 4.4 “GNSS Integrity and Quality Control”.

The precise products can be accessed at ftp://igs.gnsswhu.cn/pub/whu/phasebias/. 
Latest updates for Support, Training courses and FAQ can be found at https://pride.whu.edu.cn. 
The copyright of this package is protected by GNU General Public License (version 3). 

Relevant publications are

* Geng J, Wen Q, Zhang Q, Li G, Zhang K (2022). GNSS observable-specific phase biases for all-frequency PPP ambiguity resolution. Journal of Geodesy, 96(11):1-18. doi:10.1007/s00190-022-01602-3
* Geng J, Chen X, Pan Y, Zhao Q (2019). A modified phase clock/bias model to improve PPP ambiguity resolution at Wuhan University. Journal of Geodesy, 93(10):2053-2067. doi:10.1007/s00190-019-01301-6
* Geng J, Chen X, Pan Y, Mao S, Li C, Zhou J, Zhang K (2019). PRIDE PPP‑AR: an open‑source software for GPS PPP ambiguity resolution. GPS Solutions, 23(91):1-10. doi:10.1007/s10291-019-0888-1
* Geng J, Yang S, Guo J (2021). Assessing IGS GPS/Galileo/BDS-2/BDS-3 phase bias products with PRIDE PPP-AR. Satellite Navigation, 2(1):1-15. doi:10.1186/s43020-021-00049-9
* Geng J, Mao S. Massive GNSS network analysis without baselines: Undifferenced ambiguity resolution. J. Geophys. Res. 2021, 126(10), e2020JB021558. doi:10.1029/2020JB021558

PRIDE PPP-AR ver. 2.2 is available for:

1) Support for GPS, GLONASS, Galileo, BDS-2/3 and QZSS;
2) High-rate GNSS data processing at rates up to 50 Hz;
3) Vienna Mapping Function 1/3 (VMF1/VMF3) for troposphere modeling;
4) Second-order ionospheric correction;
5) High-dynamic mobile platforms applicable for aerial photogrammetry, ship-borne gravimetry, etc;
6) Mitigation of receiver clock jumps;
7) Processing of multi-day GNSS data, supporting continuous observations up to 32 days;
8) Satellite attitude quaternions;
9) A GUI lite-version specifically designed for early career researchers;
10) Ambiguity-float PPP with support for data dating back to 1994 when selective availablity (SA) was on;
11) GPS/Galileo/BDS-2/3 PPP-AR with Bias-SINEX format phase biases (ftp://igs.gnsswhu.cn);
12) Adoption of both rapid products (RAP) and real-time archived products (RTS) for more timeliness.

The improvements made in PRIDE PPP-AR version 2.2 include:

1. New batch script name “*pdp3*”, updated command-line parameters;
2. Multi-day processing capability up to 32 days;
3. Quaternion products support;
4. No longer requires DCB products;
5. New default products after 2020: Multi-GNSS satellite orbit, clock, bias, quaternion, and ERP products from Wuhan University;
6. “leap.sec” now required, “glonass\_chn” replaced by “sat\_parameters”;
7. GUI version with additional plotting functions;
8. Known bugs have been fixed.

`Notes: The multi-GNSS satellite orbit/clock/bias/quaternion products from Wuhan University are required by PRIDE PPP-AR version 2.2`

## Version History

### v1.0

2019-03-21

Release of PRIDE PPP-AR v1.0

### v1.1

2019-04-03

* Small bug fixing
* RINEX-3 support
* Fixed bug for high-rate computation
* Support Linux-32 system (src/lib/shard/linux-32)
* Support Mac OS system (src/lib/shard/mac)

### v1.2

2019-05-01

* Support VMF1

### v1.3

2019-05-23

* Auto-selection of IGS ATX
* Change SP3 from COD to WHU since 2019

2019-06-01

* Add src/utils/xyz2enu

2019-07-12

* Support rapid phasebias product

### v1.4

2019-07-16

* Add function: receiver clock jump check & recover
* Print table valid time by pride\_pppar
* Compatibility fixing for pride\_pppar.sh
* If 'rnx2rtkp' doesn't work, please download the source code through
   'https://github.com/tomojitakasu/RTKLIB/tree/rtklib\_2.4.3' and compile it by yourself.
   The binary we provided is a 32-bit version.

2019-09-05

* pride\_pppar.sh: small bugs fixed
* table: igs14.atx updated

2019-12-15

* install.sh: add install tips for src/lib/libpridepppar.so
* pride\_pppar.sh: fix known bugs & add error replay for debug
* table: jpleph\_de405 updated (valid until 2040-007)
* table: update IGS14.atx (igs14\_2082.atx)

### v2.0

2021-05-21

Release of PRIDE PPP-AR v2.0

### v2.1

2021-09-06

Release of PRIDE PPP-AR v2.1
 
* Support for quaternion products

### v2.2

2022-04-07

Release of PRIDE PPP-AR v2.2

* New batch script name “*pdp3*”, updated command-line parameters
* Multi-day processing capability up to 5 days
* No longer requires DCB products
* New default products after 2020: Multi-GNSS satellite orbit, clock, bias, quaternion, and ERP products from Wuhan University (WUM0MGXRAP)
* “leap.sec” now required, “glonass\_chn” replaced by “sat\_parameters”
* GUI version with additional plotting functions
* Known bugs have been fixed

2022-05-07

* install.sh: default table directory can be set outside /home
* pdp3.sh: adjusted PCO/PCV models for CODE products
* pdp3.sh: improved compatibility of RINEX observation file naming
* pdp3.sh: added support for file paths with spaces
* pdp3.sh: added support for older versions of *wget*
* table: added M14.ATX

2022-06-20

* pdp3.sh: added alerts for improper installation operations
* pdp3.sh: added OFFLINE mode to save time from *wget* calls
* pdp3.sh: added SA mode for SPP processing
* pdp3.sh: bug fixes for ANTEX file damage in batch processing
* pdp3.sh: bug fixes for unexpected break off in multi-day processing
* pdp3.sh: increased maximum number of days for multi-day processing to 32

2022-10-28

* pdp3.sh: increased maximum processing interval to 300.0 secs
* pdp3.sh: fixed syntax errors in output
* spp: aligned initial coordinates file timestamp with observations
* lsq: fixed potential fatal issues caused by rounding errors

2022-11-21

* New default products for dates between 1995 and 2020: IGS Repro3 combination products (IGS2R03FIN) that supprot PPP-AR after 2020
* pdp3.sh: increased decimal places in time range from two to three
* table: added igsR3\_2135.atx, removed M14.ATX

2023-02-14

* Added IGS20 compatilibility with new ANTEX file, SINEX file, and "APC\_MODEL" keyword
* pdp3.sh: decreased the minimum processing time range for loose editing to 30.0 secs
* pdp3.sh: priori station coordinates can now be specified in the sit.xyz manually
* lib: adopted complete phase center correction on the Melbourn-Wubbena combination
* lib: fixed the issue of unable to read BDS B1I observations from RINEX 3.02 format
* lib: fixed the issue of unable to read zero biases from bias product
* lsq: untrusted coordinates (that previously marked with “\*”) are no longer exproted to the kin file
* table: added igs20\_2247.atx

2023-02-27

* pdp3.sh: consider specifying the language/region setting to "en_US.UTF-8" to avoid issues with decimal comma
    - please contact us if further issues persist in your region
* table/config\_template: BDS GEO satellites are disabled by default to enhance processing quality

2023-03-31

* lib: fixed the issue of reading multi-day quaternions
* lib: improved the substitution rules for vacant multi-GNSS PCOs/PCVs and fixed related read issues
* redig: imporved the ability to detect minor cycle slip 
* spp: fixed the issue of unable to process BDS-only observations
 
## Getting in Touch

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FPrideLab%2FPRIDE-PPPAR&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=PAGE+VIEWS&edge_flat=false)](https://hits.seeyoufarm.com)

* You can contact us for **bug reports** and **comments** by sending an email or leaving a message on our website:
  * Email: <pride@whu.edu.cn>
  * Website: <http://pride.whu.edu.cn>
* For Chinese users, we provide Tencent **QQ Group** service.
  * QQ group: 971523302

## License

***Copyright (C) 2022 by Wuhan University, All rights reserved.***

---
