OCoCNANS
========

OCoCNANS (Outer Crust of Cold NonAccreting Neutron Stars) is a Python module 
that can be used to evaluate the ground state of matter in the outer crust of 
cold isolated neutron stars for a given a nuclear mass table.

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
and `ococnans.outer_crust`, and are used like this:

``` py 
import ococnans as oc
hfb26 = oc.read_masstable("masstables/hfb26.data", sep=' ', 
        mexcess=True, useexpdata=True)
ocrust = oc.outer_crust(hfb26, pressure_step=0.003)
```

`ococnans.outer_crust` returns a NumPy array for the ground state of matter in 
the outer crust. Columns are baryon number density in fm$^{-3}$, pressure in 
MeV, neutron chemical potential in MeV, proton chemical potential in MeV, mass 
number of equilibrium nucleus, and charge number number of equilibrium nucleus.

The NumPy array can eventually be converted to a pandas DataFrame and 
Matplotlib can be used for visualization.

``` py 
import pandas as pd
ocrust = pd.DataFrame(ocrust, columns=["nb", "pres", "mun", "mup", "aa", "zz"])

import matploblib.pyplot as plt
import matplotlib.gridspec as gridspec
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(2, 2)

ax = fig.add_subplot(gs[0, :])
ax.plot(ocrust.pres, ocrust.aa, label="$A$")
ax.plot(ocrust.pres, ocrust.zz, label="$Z$")
ax.set_xlabel("pressure [MeV/fm$^3$]")
ax.set_xlim(min(ocrust.pres), max(ocrust.pres))
ax.set_xscale("log")
ax.legend(loc='best')

ax = fig.add_subplot(gs[1, 0])
ax.plot(ocrust.pres, ocrust.mun, label="$\mu_n$")
ax.plot(ocrust.pres, ocrust.mup, label="$\mu_p$")
ax.set_xlabel("pressure [MeV/fm$^3$]")
ax.set_ylabel("chemical potential [MeV]")
ax.set_xlim(min(ocrust.pres), max(ocrust.pres))
ax.set_xscale("log")
ax.legend(loc='best')

ax = fig.add_subplot(gs[1, 1])
ax.plot(ocrust.nb, ocrust.pres)
ax.set_xlabel("baryon number density [fm$^{-3}$]")
ax.set_ylabel("pressure [MeV/fm$^3$]")
ax.set_xlim(min(ocrust.nb), max(ocrust.nb))
ax.set_xscale("log")
ax.set_yscale("log")
```

![OuterCrust_HFB-26](example.png "Ground state of matter in the outer crust for HFB-26")
