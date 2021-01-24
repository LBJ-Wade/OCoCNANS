import numpy as np
from scipy.optimize import fsolve


"""Physical constants (https://physics.nist.gov/cuu/Constants/index.html)"""
HBARC = 197.3269788 # in MeV.fm
ALPHAFS = 7.2973525664e-3 # fine-structure constant
MEL = 0.5109989461 # electron mass energy equivalent in MeV
RMN = 939.5654133 # neutron mass energy equivalent in MeV
RMP = 938.2720813 # proton mass energy equivalent in MeV
AMU = 931.4940954 # atomic mass constant energy equivalent in MeV

np.seterr(invalid="ignore")


class ElectronGas:
    """Class for the relativistic electron gas. Expressions are taken from
    [Haensel P., Potekhin A. Y., Yakovlev D. G., 2007, Neutron Stars 1: 
    Equation of state and structure. Springer, Berlin].

    Parameters
    ----------
    electron_number_density : float or array-like of floats
        Electron number density (unit is /fm^3).

    """
    def __init__(self, electron_number_density):
        self.electron_number_density = electron_number_density
        self.fermi_wave_number = (3.0 * np.pi * np.pi
                * self.electron_number_density) ** (1.0 / 3.0)
        self.xr = HBARC * self.fermi_wave_number / MEL # relativity parameter
        self.gammar = (self.xr * self.xr + 1.0) ** 0.5

    def energy_density(self):
        """Return the electron gas energy density (unit is MeV/fm^3). 
        Rest mass is included."""
        return (MEL ** 4) / 8.0 / np.pi / np.pi / (HBARC * HBARC * HBARC)\
                * ((2.0 * self.xr * self.xr + 1.0)\
                * self.xr * self.gammar - np.log(self.xr + self.gammar))

    def chemical_potential(self):
        """Return the electron chemical potential (unit is MeV). Rest mass is 
        included)."""
        return MEL * MEL * MEL / 8.0\
                / (self.fermi_wave_number * self.fermi_wave_number)\
                / (HBARC * HBARC)\
                * (self.gammar * (1.0 + 6.0 * self.xr * self.xr)\
                + self.xr *self.xr * (2.0 * self.xr * self.xr + 1.0)\
                / self.gammar - 1.0 / self.gammar)

    def pressure(self):
        """Return the electron gas pressure (unit is MeV/fm^3)."""
        return self.electron_number_density * self.chemical_potential()\
                - self.energy_density()


def lattice_energy_density(charge_number, electron_number_density):
    """Calculate the lattice energy density for a bcc structure. Expression is
    taken from [Pearson J. M., Goriely S., Chamel N., 2011, Phys. Rev. C,
    83,065810].

    Parameters
    ----------
    charge_number : int
        Atomic number of the nucleus.

    electron_number_density : float or array-like of floats
        Electron number density (unit is /fm^3).

    Returns
    -------
    eldens : float or array-like of floats
        Lattice energy density (unit is MeV/fm^3).

    """
    eldens = -0.895929255682 * (4.0 * np.pi / 3.0) ** (1.0 / 3.0)\
            * charge_number ** (2.0 / 3.0)\
            * electron_number_density ** (4.0 / 3.0)
    return eldens


def lattice_pressure(charge_number, electron_number_density):
    """Calculate the lattice pressure.

    Parameters
    ----------
    charge_number : int
        Atomic number (Z) of the nucleus.

    electron_number_density : float  or array-like of floats
        Electron number density (unit is /fm^3).

    Returns
    -------
    pl : float or array-like of floats
        Lattice pressure (unit is MeV/fm^3).

    """
    pl = lattice_energy_density(charge_number, electron_number_density) / 3.0
    return pl


def _f_pressure(x, *args):
    # args[0] (args[1]) is the charge number (pressure)
    return ElectronGas(x).pressure() + lattice_pressure(args[0], x) - args[1]


def _solve_electron_number_density(charge_number, pressure):
    electron_number_density = fsolve(_f_pressure, 
            x0=len(charge_number)*[1.0e-10],
            args=(charge_number, pressure))
    return electron_number_density


def electron_binding_energy(charge_number):
    """Return the electron binding energy for a given number of protons (unit
    is MeV). Expression is taken from [Lunney D., Pearson J. M., Thibault C., 
    2003, Rev. Mod. Phys.,75, 1021]."""
    return 1.44381e-5 * charge_number ** 2.39\
            + 1.55468e-12 * charge_number ** 5.35


def nuclear_mass_from_mass_excess(mass_number, charge_number, mass_excess):
    """Calculate the nuclear mass from the mass excess for a given nucleus.

    Parameters
    ----------
    mass number : int 
        Mass number (A) of the nucleus.

    charge_number : int
        Atomic number (Z) of the nucleus.

    mass_excess : float or array-like of floats
        Mass excess of the nucleus (unit is MeV).

    Returns
    -------
    nuclear_mass : float or array-like of floats
        Nuclear mass corresponding to the input nucleus (unit is MeV).

    """
    nuclear_mass = mass_excess + mass_number * AMU\
            - charge_number * MEL + electron_binding_energy(charge_number)
    return nuclear_mass


def gibbs_free_energy_per_nucleon(mass_number, charge_number,
        electron_number_density, nuclear_mass):
    """Calculate the Gibbs free energy per nucleon for a given electron number
    density and nucleus.

    Parameters
    ----------
    mass number : int 
        Mass number (A) of the nucleus.

    charge_number : int
        Atomic number (Z) of the nucleus.

    electron_number_density : float or array-like of floats
        Electron number density (unit is /fm^3).

    nuclear_mass : float or array-like of floats
        Nuclear mass (unit is MeV).

    Returns
    -------
    gibbs_free_en_per_nuc : float or array-like of floats
        Gibbs free energy per nucleon (unit is MeV).

    """
    wscell_volume = charge_number / electron_number_density
    baryon_number_density = mass_number / wscell_volume
    gibbs_free_en_per_nuc = nuclear_mass / mass_number\
            + 4.0 / 3.0 * lattice_energy_density(charge_number,
                    electron_number_density) * wscell_volume / mass_number\
            + charge_number / mass_number\
            * (ElectronGas(electron_number_density).chemical_potential() - MEL)
    return gibbs_free_en_per_nuc


def read_masstable(filepath, mexcess=True, sep=' ', useexpdata=False):
    """Read a nuclear mass table to be used in the calculation.

    Parameters
    ----------
    filepath : str
        Path of nuclear mass table. Columns must be proton number, neutron
        number, and either nuclear mass or mass excess in MeV.

    mexcess : bool, default=True
        If True, convert mass excess into nuclear mass.

    sep : str, default=' '
        String used to separate values. 

    useexpdata : bool, default=False
        If True, the latest experimental masses will be used when available 
        (Atomic Mass Evaluation 2016 and masses of neutron-rich copper isotopes
        from Phys. Rev. Lett. 119, 192502).

    Returns
    -------
    table : NumPy array
        Nuclear mass table. Columns are proton number, neutron number, and
        nuclear mass (in MeV).

    """
    masstable = np.loadtxt(filepath, delimiter=sep)
    if useexpdata == True:
        try: 
            masstable_with_expdata = np.loadtxt(
                    "masstables/ame2016+welker2017.data")
        except OSError:
            print(
                    "Cannot find masstables/ame2016+welker2017.data. "
                    "Experimental data will not be used.\n"
                    "Please download the file at "
                    "https://raw.githubusercontent.com/thomascarreau/OCoCNANS/"
                    "master/masstables/ame2016+welker2017.data."
                    )
        else:
            ame_nuclei = list(zip(masstable_with_expdata[:, 0],
                masstable_with_expdata[:, 1]))
            for nucleus in masstable:
                if (nucleus[0], nucleus[1]) not in ame_nuclei:
                    masstable_with_expdata = np.vstack(
                            [masstable_with_expdata, 
                                np.array(
                                    [nucleus[0], nucleus[1], nucleus[2]])])
            masstable = masstable_with_expdata[
                    (masstable_with_expdata[:, 0] > 0) & 
                    (masstable_with_expdata[:, 1] > 0)]
    if mexcess == True:
        masstable[:, 2] = nuclear_mass_from_mass_excess(masstable[:, 0] +
                masstable[:, 1], masstable[:, 0], masstable[:, 2])
    return masstable


def outer_crust_composition_at_fixed_pressure(pressure, masstable):
    """Calculate the outer crust composition at fixed pressure for a given mass
    table.

    Parameters
    ----------
    pressure : float
        Pressure (unit is MeV/fm^3).

    masstable : NumPy array
        Nuclear mass table.

    Returns
    -------
    nb_eq : float
        Equilibrium value of baryon number density (unit is /fm^3).

    gmin : float
        Equilibrium value of Gibbs free energy per nucleon (unit is MeV).

    aa_eq : int
        Equilibrium value of mass number.

    zz_eq : int
        Equilibrium value of charge number.

    """
    # calculate the electron density for atomic numbers in mass table
    zz = np.arange(np.min(masstable[:, 0]), np.max(masstable[:, 0]) + 1, 1)
    ne = _solve_electron_number_density(zz, pressure)
    ne_table = []
    for zzi in masstable[:, 0]:
        ne_table.append(ne[np.where(zz == zzi)[0][0]])
    # calculate the gibbs free energy per nucleon for each nucleus
    aa = masstable[:, 0] + masstable[:, 1]
    gibbs_free_en_per_nuc = gibbs_free_energy_per_nucleon(
            aa, masstable[:, 0], np.array(ne_table), masstable[:, 2])
    # determine ground state
    idx_sol = np.where(gibbs_free_en_per_nuc ==
            np.min(gibbs_free_en_per_nuc))[0][0] # equilibrium solution index
    aa_eq = aa[idx_sol]
    zz_eq = masstable[:, 0][idx_sol]
    nb_eq = aa_eq / zz_eq * ne_table[idx_sol]
    return nb_eq, gibbs_free_en_per_nuc[idx_sol], aa_eq, zz_eq


def outer_crust(masstable, pressure_init=9.0e-12, pressure_step=0.05):
    """Calculate the ground state of matter in the outer crust of
    cold-catalyzed neutron stars for a given mass model.

    Parameters
    ----------
    masstable : NumPy array
        Nuclear mass table.

    pressure_init : float, default=9.0e-12
        Initial value of pressure (unit is MeV/fm^3).

    pressure_step : float, default=0.05
        Normalized step of pressure (unit is MeV/fm^3). At each iteration, 
        pressure += pressure_step * pressure. While 0.05 may be sufficient to 
        have a good overview of the composition, it is recommended to set a 
        lower pressure step (< 0.005) to investigate the presence of thin 
        layers and get a proper estimation of the neutron-drip point.

    Returns
    -------
    ocrust : NumPy array
        Ground state of matter in the outer crust. Columns are baryon number
        density (in /fm^3), pressure (in MeV/fm^3), neutron chemical potential
        (in MeV), proton chemical potential (in MeV), mass number of 
        equilibrium nucleus, and charge number of equilibrium nucleus.

    Examples
    --------
    >>> import ococnans as oc
    >>> hfb26 = oc.read_masstable("masstables/hfb26.data", sep=' ',
    ...     mexcess=True, useexpdata=True)
    >>> ocrust = oc.outer_crust(hfb26, pressure_step=0.003)
    >>> import pandas as pd
    >>> ocrust = pd.DataFrame(ocrust, 
    ...     columns=["nb", "pres", "mun", "mup", "aa", "zz"])

    """
    pressure = pressure_init
    gibbs_free_en_per_nuc = -1.0e9
    ocrust = []
    while gibbs_free_en_per_nuc < RMN: # condition for neutron drip
        # calcuate equilibrium configuration at fixed pressure
        comp = outer_crust_composition_at_fixed_pressure(pressure, masstable)
        baryon_number_density = comp[0]
        gibbs_free_en_per_nuc = comp[1]
        aa = comp[2]
        zz = comp[3]
        # calculate proton chemical potential
        egas = ElectronGas(zz / aa * baryon_number_density)
        mup = gibbs_free_en_per_nuc - (egas.chemical_potential() - MEL)
        # append composition into the array
        ocrust.append([baryon_number_density, pressure, gibbs_free_en_per_nuc, 
            mup, int(aa), int(zz)])
        # iterate pressure
        pressure += pressure * pressure_step
    return np.array(ocrust)
