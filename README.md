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

PRIDE PPP-AR ver. 2.2 are available for:  

1) GPS, GLONASS, Galileo, BDS-2/3 and QZSS capable;  
2) High-rate GNSS data processing of up to 50 Hz;
3) Vienna Mapping Function 1/3 (VMF3) for troposphere modeling;  
4) Second-order ionospheric correction;  
5) High-dynamic mobile platforms applicable for aerial photogrammetry, ship-borne gravimetry, etc;  
6) Receiver clock jump mitigation;  
7) Multi-day processing;  
8) Satellite attitude quaternions;  
9) A GUI lite-version provided for very early career researchers;  
10) Ambiguity-float PPP using data dating back to 1994 when SA was on;  
11) GPS/Galileo/BDS-2/3 PPP-AR in the case of the bias-SINEX format phase biases (ftp://igs.gnsswhu.cn);
12) Adopt both rapid products (RAP) and real-time products (RTS) for more timeliness.

The modifications leading to Version 2.2 include:

1. Batch script name changed from “*pride_pppar*” to “*pdp3*”, corresponding command line input parameters also changed;
2. Support multi-day processing;
3. Support for quaternion products;
4. No more DCB products required;
5. The default products after 2020 changed to the multi-GNSS satellite orbit, clock, bias, quaternion and ERP products, which are computed and released by Wuhan University;
6. The table file “leap.sec” needs to be downloaded now, and the “glonass\_chn” table file is removed and replaced by the “sat\_parameters” table file;
7. GUI version of PRIDE PPP-AR with additional plotting functions;
8. Fix known bugs.

`Notes: The multi-GNSS multi-GNSS satellite orbit/clock/bias/quaternion products, which are computed and released by Wuhan University, are required by PRIDE PPP-AR ver. 2.2`

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
* Print table valid time by pride_pppar
* Compatibility fixing for pride_pppar.sh
* If 'rnx2rtkp' doesn't work, please download the source code through
   'https://github.com/tomojitakasu/RTKLIB/tree/rtklib_2.4.3' and compile it by yourself.
   The binary we provided is a 32-bit version.

2019-09-05

* pride_pppar.sh: small bugs fixed
* table: igs14.atx updated

2019-12-15

* install.sh: add install tips for src/lib/libpridepppar.so
* pride_pppar.sh: fix known bugs & add error replay for debug
* table: jpleph_de405 updated (valid until 2040-007)
* table: update IGS14.atx (igs14_2082.atx)

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
 
* Batch script name changed from “*pride_pppar*” to “*pdp3*”, corresponding command line input parameters also changed
* Support multi-day processing
* No more DCB products required
* The default products after 2020 changed to the multi-GNSS satellite orbit, clock, bias, quaternion and ERP products, which are computed and released by Wuhan University
* The table file “leap.sec” needs to be downloaded now, and the “glonass_chn” table file is removed and replaced by the “sat_parameters” table filex
* GUI version of PRIDE PPP-AR with additional plotting functions
* Fix known bugs

2022-05-07

* install.sh: default table directory can be set outside /home directory
* pdp3.sh: add notices for some improper operations
* pdp3.sh: add "OFFLINE" mode
* pdp3.sh: improve recognition ability of RINEX file name in multi-day processing
* pdp3.sh: support file path with spaces
* pdp3.sh: support old versions of wget
* pdp3.sh: bug fixes for not using M14.ATX in multi-day processing with CODE products, adjust PCO/PCV models used for different CODE products
* table: add M14.ATX
 
## Getting in Touch

* You can contact us for **bug reports** and **comments** by sending an email or leaving a message on our website:
  * Email: <pride@whu.edu.cn>
  * Website: <http://pride.whu.edu.cn>
* For Chinese users, we provide Tencent **QQ Group** service.
  * QQ group: 971523302

## License

***Copyright (C) 2022 by Wuhan University, All rights reserved.***

---
