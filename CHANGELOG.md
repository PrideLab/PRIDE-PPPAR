# Change Log

## PRIDE PPP-AR version 3.2.0

### 2025-10-20 (v3.2.0)
* `arsig`: A machine learning–based model is adopted by default for ambiguity validation
* `tedit`: Single constellation data edition is allowed
* `pdp3`: A new “-wcc” option has been added to the command line to support the use of Wuhan Combination Center combined products
* `pdp3`: The default FTP address for WUM0MGXRAP products is changed to "ftps://bdspride.com/wum/"

### 2025-09-06 (v3.1.5)
* `lsq`: add 'isb' estimation for receiver clocks  

### 2025-05-17 (v3.1.4)
* `gui`: update rinex2 reading function and high-order ionosphere correction function

### 2025-04-11 (v3.1.3)
* `lib` & `orbit` : optimize reading RINEX-2.12 data and orbit product code, respectively

### 2025-02-25 (v3.1.2)
* `lib`: update IGRF12.F

### 2025-01-25 (v3.1.2)
* `lsq`: enable the processing of RINEX files with non-standard time label

### 2025-01-19 (v3.1.1)
* `pdp3`: Disable the use of RTS products in MHM modeling

### 2025-01-06 (v3.1)
* `mhm`: a new module enabling multipath delay compensation based on Multipath Hemispherical Map model (MHM)
* `gui`: add support for RINEX versions greater than 4.0 


### 2024-09-08 (v3.0)
* `lib` & `tedit`: use DOCB block in BIAS-SINEX file to decide whether reset ambiguity at midnight
* `otl`: enable self-computed ocean tide loading correction in case of no oceanload coefficients for user stations
* `lsq`: add a function to support random walk constraint between epochs in position domain (P mode)

### 2024-04-08 (v3.0)

* `lib` & `spp`: enable multi-day processing of non-consecutive daily RINEX files
* `tedit`: improve implementation of loose editing for kinematic positioning
* `lsq`: reduce the weight of GLONASS pseudorange observations to improve multi-GNSS positioning performance

### 2024-02-08 (v3.0)

* Employ the latest **WUM0MGXRTS** products to perform all-frequency PPP-AR on GNSS data with hours latency
  * Real-time archived from Wuhan University GPS/Galileo/BDS phase bias stream and updated every 3 hours
* `lib` & `spp`: fix bug for unexpected termination in multi-day processing

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

1. `pdp3_mac.sh` has been integrated into `pdp3.sh`
2. Older versions of PRIDE PPP-AR cannot recognize the new namings of WUM0MGXRAP and WUM0MGXRTS products
3. Older versions of PRIDE PPP-AR cannot read the new settings in configuration file
4. Older versions of PRIDE PPP-AR cannot read the new format in result files

## PRIDE PPP-AR version 2

### 2023-09-28 (v2.2)

* `pdp3`: specify the language/region setting to `en_US.UTF-8` to prevent problems with decimal comma
  * please contact us if you still encounter problems in your region
* `pdp3`: support new GIM products for 2nd ionospheric corrections from 2022/330
* `pdp3`: fix bug for not recognizing some RINEX observation files
* `lib`: fix bug for reading multi-day quaternions
* `lib`: improve the substitution rules for missing multi-GNSS PCOs/PCVs
* `spp` & `tedit`: fix bug for not processing BDS-only observations
* `table`: disable BDS GEO satellites in `config_template` by default to enhance processing quality

### 2023-02-14 (v2.2)

* Add **IGS20** compatilibility with new ANTEX file, SINEX file, and `APC_MODEL` keyword
* `pdp3`: reduce the minimum processing time range for loose editing to 30.0 seconds
* `pdp3`: allow manual specification of priori station coordinates in `sit.xyz` in `F` mode
* `lib`: apply complete APC corrections on the Melbourn-Wubbena combination instead of approximate ones
* `lib`: fix bug for not reading BDS B1I observations from RINEX 3.02 format
* `lib`: fix bug for not reading zero biases from code/phase OSB product
* `table`: add `igs20_2247.atx` for products in IGS20 reference frame

### 2022-11-21 (v2.2)

* Change default products from CODE to **IGS2R03FIN** for dates between 1995 and 2020
  * IGS Repro3 combination products that supprot PPP from 1995 and PPP-AR from 2000
* `pdp3`: increase the number of decimal digits in time range to 3
* `table`: add `igsR3_2135.atx` for Repro3 products, remove `M14.ATX` for CODE products

### 2022-10-28 (v2.2)

* `pdp3`: extend maximum processing interval to 300.0 seconds
* `pdp3`: correct syntax errors in output
* `spp`: match the timestamp of initial coordinates with observations
* `lsq`: fix potential fatal issues caused by rounding errors

### 2022-06-20 (v2.2)
* `install.sh`: add warnings for improper installation operations
* `pdp3`: add OFFLINE mode to avoid calling `wget` in an off-line computer
* `pdp3`: fix bug for corrupting IGS ANTEX file in batch processing
* `pdp3`: fix bug for unexpected termination in multi-day processing
* `pdp3`: enhance capability of **multi-day processing** to 32 days
* `spp`: support SA mode for dates before 2000

### 2022-05-07 (v2.2)
* `install.sh`: allow default table directory to be set outside `/home`
* `pdp3`: enhance compatibility of RINEX observation file naming
* `pdp3`: support file paths with spaces
* `pdp3`: support older versions of `wget` (without `--show-progress`)
* `table`: add `M14.ATX` for CODE products

### 2022-04-07 (v2.2)

Release of **PRIDE PPP-AR v2.2**

* Change command-line script name from `pride_pppar` to `pdp3`
* Redesign interactions with command-line interface without fixed arguments
* Employ **WUM0MGXRAP** products from Wuhan University for dates from 2020
  * Rapid multi-GNSS satellite orbit, clock, bias, quaternion, and ERP products from Wuhan University
  * Apply code/phase OSB products in Bias-RINEX format instead of DCB products
* Employ ***WUM0MGXRTS** products from Wuhan University for ultra-rapid processing
* Support continuous **multi-day processing** up to 5 days
* Release a GUI version for Windows/MacOS with interactive plotting functions
* Add `src/utils/pbopos` to convert `pos` file to PBO GPS Time Series Format
* `table`: add `config_template` to save command-line settings
* `table`: add `leap.sec` to synchronize leap seconds from IERS
* `table`: replace `glonass_chn` with `sat_parameters`

### 2021-09-06 (v2.1)

Release of **PRIDE PPP-AR v2.1**

* Support satellite attitude quaternion products

### 2021-05-21 (v2.0)

Release of **PRIDE PPP-AR v2.0**

* Support **Multi-GNSS:** GPS, GLONASS, Galileo, BDS-2/3, QZSS

## PRIDE PPP-AR version 1

### 2019-12-15 (v1.4)

* `install.sh`: add install tips for `src/lib/libpridepppar.so`
* `pride_pppar`: fix known bugs & add error replay for debug
* `table`: update `jpleph_de405` (valid until 2040-007)

### 2019-09-05 (v1.4)

* `pride_pppar`: fix small bugs

### 2019-07-16 (v1.4)

* Add function: receiver clock jump check & recover
* Print table valid time by `pride_pppar`
* Enhance compatibility of `pride_pppar`

### 2019-07-12 (v1.3)

* Support rapid phasebias product

### 2019-06-01 (v1.3)

* Add `src/utils/xyz2enu`

### 2019-05-23 (v1.3)

* Enable auto-selection of IGS ANTEX file
* Change defult SP3 products from COD to WHU for the dates since 2019

### 2019-05-01 (v1.2)

* Support VMF1

### 2019-04-03 (v1.1)

* Fix small bugs
* Fix bug for high-rate computation
* Support RINEX-3 
* Support Linux-32 system (`src/lib/shard/linux-32`)
* Support Mac OS system (`src/lib/shard/mac`)

### 2019-03-21 (v1.0)

Release of **PRIDE PPP-AR v1.0**
