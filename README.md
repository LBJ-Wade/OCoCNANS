OCoCNANS
========

OCoCNANS, for Outer Crust of Cold Non-Accreting Neutron Stars, is a small 
Python code for calculating the outer crust composition and energetics of a 
cold non-accreting neutron star for a given nuclear mass table.

Requirements
------------

* Python 3

* [NumPy](https://numpy.org/install/)

* [SciPy](https://scipy.org/install.html)

Getting started
---------------

    git clone https://github.com/thomascarreau/OCoCNANS
    cd OCoCNANS

The only two user-facing functions in the module are `ococnans.read_masstable` 
and `ococnans.outer_crust`:

``` { .py }
import ococnans as oc
hfb26 = oc.read_masstable("masstables/hfb26/data", sep=' ', 
    mexcess=True, useexpdata=True)
ocrust = oc.outer_crust(hfb26, pressure_step=0.003)
```

> ### Available mass tables

> The format of the mass tables is Z, N, mass excess (`deps`).

> * AME1995 (with and without interpolated values)

> * AME2003 (with and without interpolated values)

> * AME2012 (with and without interpolated values)

> * BR2013

> * ETFSI12

> * FRDM95

> * HFB14

> * HFB26

> * KTUY05

> * SVM13

> * TCSM12

> * TCSM13
