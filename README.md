# OCoCNANS

OCoCNANS, for Outer Crust of Cold Non-Accreting Neutron Stars, is a Python code for calculating the outer crust composition and energetics of a cold non-accreting neutron star for a given nuclear mass table.

## Requirements 

* Python 2.7

* [SciPy](https://www.scipy.org/scipylib/index.html)


## Usage

    python2.7 ococnans.py mass_table.data ocrust.out

or

    ./ococnans.py mass_table.data ocrust.out

if you set execute permission with

    chmod +x ococnans.py

### Available mass tables

* AME1995 (with and without interpolated values)

* AME2003 (with and without interpolated values)

* AME2012 (with and without interpolated values)

* BR2013

* ETFSI12

* FRDM95

* HFB14

* HFB26

* KTUY05

* SVM13

* TCSM12

* TCSM13

See the name of the associated files in `mass_tables/`.
