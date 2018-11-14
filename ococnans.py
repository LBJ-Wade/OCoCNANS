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
AMU = 931.4940954

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

def get_electron_density(zz_, pp_):
    ne = float(fsolve(lambda ne: f_pressure(zz_, ne, pp_), 1.e-10))
    return ne

def electron_binding_energy(zz_):
    return 1.44381e-5*pow(zz_,2.39) + 1.55468e-12*pow(zz_,5.35)

def nuclear_mass(aa_, zz_, deps_):
    return deps_ + aa_*AMU - zz_*MEL + electron_binding_energy(zz_)

def gibbs_free_energy_per_nucleon(aa_, zz_, ne_, deps_):
    vws = zz_/ne_
    nb = aa_/vws
    e = nuclear_mass(aa_, zz_, deps_)
    return e/aa_ + 4./3.*lattice_energy_density(zz_, ne_)*vws/aa_ \
            + zz_/aa_*electron_chemical_potential(ne_)

def get_outer_crust_composition(pp_, mass_table):
    gmin = 1.e99
    f_data = open(mass_table, 'r')
    for line in f_data:
        columns = line.strip().split()
        zz = float(columns[0])
        nn = float(columns[1])
        deps = float(columns[2]) # mass excess
        if (nn > zz and nn % 2 == 0 and zz % 2 == 0):
            ne = get_electron_density(zz, pp_)
            aa = zz + nn
            g = gibbs_free_energy_per_nucleon(aa, zz, ne, deps)
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
    print "Layer \t n_B(max) \t P(max) \t A \t Z"
    print "---------------------------------------------------"
    comp = get_outer_crust_composition(5.9e-13, mass_table)
    aa_sav = int(comp[2])
    zz_sav = int(comp[3])
    pp = 6.e-13
    nb = 0.
    nlayer = 0
    f_ocrust = open(outfile, 'w')
    while (nb < 3.e-4):
        comp = get_outer_crust_composition(pp, mass_table)
        nb = comp[0]
        pmax = pp
        nbmax = nb
        if ("ame" in mass_table and nb > 6.e-5):
            break
        g = comp[1]
        aa = int(comp[2])
        zz = int(comp[3])
        f_ocrust.write(str(nb) + "\t" + str(pp) + "\t" \
                + str(g) + "\t" + str(aa) + "\t" + str(zz) + "\n")
        if (aa != aa_sav or zz != zz_sav):
            nlayer += 1
            print str(nlayer) + " \t " + '{:0.3e}'.format(nbmax) \
                    + " \t " + '{:0.3e}'.format(pmax) + " \t " \
                    + str(aa_sav) + " \t " + str(zz_sav)
            aa_last_layer = aa_sav
            zz_last_layer = zz_sav
        aa_sav = aa
        zz_sav = zz
        pp += pp/4.
    f_ocrust.close()
    if(aa_sav != aa_last_layer or zz_sav != zz_last_layer):
        nlayer += 1
        print str(nlayer) + " \t " + "---------" \
                + " \t " + "---------" + " \t " \
                + str(aa_sav) + " \t " + str(zz_sav)
    print ""

if __name__=="__main__":
    main()
