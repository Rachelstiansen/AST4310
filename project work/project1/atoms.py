import numpy
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy import units
from astropy import constants

# Plot defaults and quantity supportS
from astropy.visualization import quantity_support
quantity_support()

class Atom:
    """
    Reads atomic data, calculates level populations according to Boltzmann's law,
    and ionisation fractions according to Saha's law.
    """
    
    def __init__(self, atomfile=None):
        """
        Parameters
        ----------
        atomfile : string, optional
            Name of file with atomic data. If not present, atomic data needs
            to be loaded with the .read_atom method.
        """
        self.loaded = False
        if atomfile:
            self.read_atom(atomfile)
        
    def read_atom(self, filename):
        """
        Reads atom structure from text file.
        
        Parameters
        ----------
        filename: string
            Name of file with atomic data.
        """
        tmp = numpy.loadtxt(filename, unpack=True)
        self.n_stages = int(tmp[2].max()) + 1
        # Get maximum number of levels in any stage
        self.max_levels = 0
        for i in range(self.n_stages):
            self.max_levels = max(self.max_levels, (tmp[2] == i).sum())
        # Populate level energies and statistical weights
        # Use a square array filled with NaNs for non-existing levels
        chi = numpy.empty((self.n_stages, self.max_levels))
        chi.fill(numpy.nan)
        self.g = numpy.copy(chi)
        for i in range(self.n_stages):
            nlevels = (tmp[2] == i).sum()
            chi[i, :nlevels] = tmp[0][tmp[2] == i]
            self.g[i, :nlevels] = tmp[1][tmp[2] == i]
        # Put units, convert from cm-1 to Joule
        chi = (chi / units.cm).to('aJ', equivalencies=units.spectral())
        # Save ionisation energies, saved as energy of first level in each stage
        self.chi_ion = chi[:, 0].copy()
        # Save level energies relative to ground level in each stage
        self.chi = chi - self.chi_ion[:, numpy.newaxis]
        self.loaded = True
        
    def compute_partition_function(self, temperature):
        """
        Computes partition functions using the atomic level energies and
        statistical weights.
        
        Parameters
        ----------
        temperature: astropy.units.quantity (scalar or array)
            Gas temperature in units of K or equivalent.
        """
        temp = temperature[numpy.newaxis, numpy.newaxis]
        result = numpy.nansum(
            self.g[..., numpy.newaxis] * 
            numpy.exp(-self.chi[..., numpy.newaxis]/constants.k_B / temp),
            axis=1
        )
        return result
    
    def compute_excitation(self, temperature):
        """
        Computes the level populations relative to the ground state,
        according to the Boltzmann law.
        
        Parameters
        ----------
        temperature: astropy.units.quantity (scalar or array)
            Gas temperature in units of K or equivalent.
        """
        pfunc = self.compute_partition_function(temperature)
        # reshape arrays to allow broadcast
        temp = temperature[numpy.newaxis, numpy.newaxis]
        g_ratio = self.g[..., numpy.newaxis] / pfunc[:, numpy.newaxis]
        chi = self.chi[..., numpy.newaxis]
        return g_ratio * numpy.exp(-chi / constants.k_B / temp)
       
    def compute_ionisation(self, temperature, electron_pressure):
        """
        Computes ionisation fractions according to the Saha law.
        
        Parameters
        ----------
        temperature: astropy.units.quantity (scalar or array)
            Gas temperature in units of K or equivalent.
        electron_pressure: astropy.units.quantity (scalar)
            Electron pressure in units of Pa or equivalent.
        """
        N_e = electron_pressure / (constants.k_B * temperature)
        pfunc = self.compute_partition_function(temperature)
        saha_const = ((2 * numpy.pi * constants.m_e * constants.k_B * temperature) / (constants.h**2))**(3/2)
        nstage = numpy.zeros_like(pfunc) / units.m**3
        nstage[0] = 1. / units.m**3
        for r in range(self.n_stages -1):
            nstage[r + 1] = (nstage[r] / N_e * 2 * saha_const * pfunc[r + 1] / pfunc[r] 
                             * numpy.exp(-self.chi_ion[r + 1, numpy.newaxis] / (constants.k_B * temperature[numpy.newaxis])))
        
        # nansum is needed because the last stage may have only one level
        # (e.g. only ionisation potential)
        return nstage / numpy.nansum(nstage, axis=0) # get number of

    def compute_populations(self, temperature, electron_pressure):
        """
        Computes relative level populations for all levels and all
        ionisation stages using the Bolzmann and Saha laws.
        
        Parameters
        ----------
        temperature: astropy.units.quantity (scalar or array)
            Gas temperature in units of K or equivalent.
        electron_pressure: astropy.units.quantity (scalar)
            Electron pressure in units of Pa or equivalent.
        """
        return (self.compute_excitation(temperature) * self.compute_ionisation(temperature, electron_pressure)[:, numpy.newaxis])

    def plot_payne(self, temperature, electron_pressure):
        """
        Plots the Payne curves for the current atom.
        
        Parameters
        ----------
        temperature: astropy.units.quantity (array)
            Gas temperature in units of K or equivalent.
        electron_pressure: astropy.units.quantity (scalar)
            Electron pressure in units of Pa or equivalent.
        """
        # calculate population (relative) for a given temperature and pressure
        pops = self.compute_populations(temperature, electron_pressure)
            
        fig, ax = plt.subplots()
        ax.plot(numpy.tile(temperature, (self.n_stages, 1)).T, pops[:, 0].T, "b-") # temp vs pop of ground state of all ionization stages
        n_levels = self.chi.shape[1] # number of levels
        
        if n_levels > 1:
            ax.plot(numpy.tile(temperature, (self.n_stages, 1)).T, pops[:, 1].T, "r--") # temp vs first ionization state, second etc.
        if n_levels > 2:
            ax.plot(numpy.tile(temperature, (self.n_stages, 1)).T, pops[:, 2].T, "k:")
        ax.set_yscale("log")
        ax.set_ylim(1e-6, 1.1)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Populations")