#!/usr/bin/python2.7

import sys
from math import pi,log,pow
from scipy.optimize import fsolve

PI2 = 2.*pi

HBARC = 197.3269788
ALPHAFS = 7.2973525664e-3
MEL = 0.5109989461
RMN = 939.5654133
RMP = 938.2720813

def electron_energy_density(ne_):
    xr = HBARC*(3.*PI2*ne_)**(1./3.)/MEL
    gammar = pow(xr*xr+1.,1./2.)
    return pow(MEL,4.)/8./PI2/pow(HBARC,3.)\
            *((2.*xr*xr+1.)*xr*gammar - log(xr + gammar)) 

def electron_chemical_potential(ne_):
    xr = HBARC*(3.*PI2*ne_)**(1./3.)/MEL
    gammar = pow(xr*xr+1.,1./2.)
    return pow(MEL,3.)/8./pow(3.*ne_*PI2,2./3.)\
            /pow(HBARC,2.)*(gammar*(1.+6.*xr*xr) \
            + xr*xr*(2.*xr*xr+1.)/gammar - 1./gammar)

def electron_pressure(ne_):
    return ne_*electron_chemical_potential(ne_) \
            - electron_energy_density(ne_)

def lattice_energy_density(zz_, ne_):
    return -0.895929255682*pow(4.*pi/3.,1./3.)\
            *pow(zz_,2./3.)*pow(ne_,4./3.)

def lattice_pressure(zz_, ne_):
    return lattice_energy_density(zz_, ne_)/3.

def f_pressure(zz_, ne_, pp_):
    return electron_pressure(ne_) + lattice_pressure(zz_, ne_) - pp_

def get_electron_density(pp_, zz_):
    ne = float(fsolve(lambda ne: f_pressure(zz_, ne, pp_), 1.e-10))
    return ne

def gibbs_free_energy_per_nucleon(aa_, zz_, ne_, b_, pp_):
    vws = zz_/ne_
    nb = aa_/vws
    return b_ + 4./3.*lattice_energy_density(zz_, ne_)*vws/aa_ \
            + zz_/aa_*electron_chemical_potential(ne_) \
            + zz_/aa_*(RMP-RMN) + RMN

def get_outer_crust_composition(pp_, mass_table):
    gmin = 1.e99
    f_data = open(mass_table, 'r')
    for line in f_data:
        columns = line.strip().split()
        zz = float(columns[0])
        aa = float(columns[1])
        b = float(columns[2])
        ne = get_electron_density(pp_, zz)
        g = gibbs_free_energy_per_nucleon(aa, zz, ne, b, pp_)
        if g < gmin:
            aa_eq = aa
            zz_eq = zz
            gmin = g
    f_data.close()
    return aa/zz*ne, gmin, aa_eq, zz_eq
    

def main():
    mass_table = sys.argv[1]
    outfile = sys.argv[2]
    print ""
    print "GROUND STATE OF MATTER IN THE OUTER CRUST AT ZERO TEMPERATURE"
    print "============================================================="
    print "Mass table: " + mass_table
    print "Output file: " + outfile + " (format is nB, P, g, A, Z)"
    print ""
    print "Layer \t n_B(min) \t P(min) \t A \t Z"
    print "---------------------------------------------------"
    pp = 6.e-13
    nb = 0.
    aa_sav = 0
    zz_sav = 0
    nlayer = 0
    pmin = 0
    nbmin = 0
    f_ocrust = open(outfile, 'w')
    while(nb < 1.e-5):
        comp = get_outer_crust_composition(pp, mass_table)
        nb = comp[0]
        g = comp[1]
        aa = int(comp[2])
        zz = int(comp[3])
        f_ocrust.write(str(nb) + "\t" + str(pp) + "\t" \
                + str(g) + "\t" + str(aa) + "\t" + str(zz) + "\n")
        if (aa != aa_sav or zz != zz_sav):
            nlayer += 1
            print str(nlayer) + " \t " + '{:0.3e}'.format(nbmin) \
                    + " \t " + '{:0.3e}'.format(pmin) + " \t " \
                    + str(aa) + " \t " + str(zz)
            nbmin = nb
            pmin = pp
        aa_sav = aa
        zz_sav = zz
        pp += pp/2.
    f_ocrust.close()
    print ""

if __name__=="__main__":
    main()
