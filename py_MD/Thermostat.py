from .Integrator import Integrator
from types import FunctionType
import numpy as np

class Thermostat(Integrator):
    """
    Abstract class from which all thermostats will inherit.
    Notice that a thermostat is an integrator of the equations-of-motion with
    special modifications to ensure sampling of the NVT ensemble.
    Args:
        geometry    : initial geometry to begin dynamics from
        masses      : masses of all atoms given in geometry
        potential   : potential function (likely derived from a potential object)
        dt          : time step for the integration (in atomic units)
        temperature : target temperature of the simulation
    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature: float):
        super().__init__(geometry, masses, potential_function, dt, temperature=temperature)

    def velocity_verlet_integrate(self):
        """
        Velocity verlet integration provided to minimize copy-and-pasting without
        having to do something like inherit directly from the Velocity-Verlet class
        which should remain logically distinct from Thermostat objects because Integrators
        belong to the NVE ensemble.
        """
        # half-step of velocities
        self.current_velocities += 0.5 * self.current_accelerations * self.dt

        # update geometry
        self.current_geometry += self.current_velocities * self.dt
        
        self.update_energy_and_accelerations()

        # second half-step of velocities
        self.current_velocities += 0.5 * self.current_accelerations * self.dt

    def integrate(self):
        raise NotImplementedError("You tried to integrate with the base class thermostat. This is an abstract base class. Please choose a child of this class instead.")

################################################################################################

class Velocity_Rescaling(Thermostat):
    """
    Implements a thermostat which simply re-scales velocities to match the desired temperature.
    Note that this is known to NOT sample the canonical ensemble, and thus should be avoided for
    any real NVT simulations.
    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature: float):
        super().__init__(geometry, masses, potential_function, dt, temperature)
    
    def integrate(self):
        """
        Removes center-of-mass motion from system and then re-scales velocities.
        Propagates using velocity-verlet integration.
        """
        self.velocity_verlet_integrate()
        self.rescale_velocities()
        self.remove_com_motion()

################################################################################################

class Andersen_Thermostat(Thermostat):
    """
    Implements the Andersen theromstat which works by modelling collisions with a heat bath as
    a Poisson process. When these stochastic collisions occur, the frequency of which is governed
    by the parameter nu, the velocities are resampled from Maxwell-Boltzmann distribution.
    Because the velocities are discontinuous, no dynamical information can be extracted from trajectories
    which use this thermostat.
    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature: float,
                nu=0.1):
        super().__init__(geometry, masses, potential_function, dt, temperature)
        self.nu = nu
    
    def integrate(self):
        self.velocity_verlet_integrate()
        self.resample_velocities()
        self.center_positions()

    def resample_velocities(self):
        """
        Draw uniform random numbers to determine if a collision occurs. If collision occurs,
        then resample velocity from maxwell-boltzmann.
        """
        velocities_to_resample = list(*np.nonzero(self.dt * self.nu > np.random.uniform(size=len(self.current_velocities))))
        new_velocities = self.sample_maxwell_boltzmann_velocities(velocities_to_resample)
        for i, i_velocity in enumerate(velocities_to_resample):
            self.current_velocities[i_velocity, :] = new_velocities[i, :]

################################################################################################

class Langevin_Thermostat_GJF(Thermostat):
    """
    Implements the Langevin theromstat which models the system via Langevin dynamics. This is
    a stochastic thermostat, so is not time-reversible, but correctly samples the canonical 
    ensemble and also provides valid dynamics. This particular Langevin thermostat is known
    to underestimate the kinetic temperature. That is, the time-dynamical sampling of kinetic
    properties is not strictly correct, but the configurational one is.

    This is adapted from:
    Farago, "A simple and effective Verlet-type algorithm for simulating Langevin dynamics"
    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature: float,
                alpha=25.0,
                mass_weighted=True):
        super().__init__(geometry, masses, potential_function, dt, temperature)
        self.mass_weighted = mass_weighted
        if self.mass_weighted:
            self.alpha = alpha * self.masses
        else:
            self.alpha = alpha

        self.a = (1 - self.alpha * self.dt / (2 * self.masses)) / (1 + self.alpha * self.dt / (2 * self.masses))
        self.b = 1 / (1 + self.alpha * self.dt / (2 * self.masses))
        self.sigma = np.sqrt(2 * self.alpha * self.temperature * self.dt)
    
    def integrate(self):
        # get the random numbers and copies of forces from previous step
        if self.mass_weighted:
            gaussian_random_numbers = np.zeros(self.current_geometry.shape)
            for i in range(len(gaussian_random_numbers)):
                gaussian_random_numbers[i,:] = np.random.default_rng().normal(0.0, self.sigma[i], size=3)
        else:
            gaussian_random_numbers = np.random.default_rng().normal(0.0, self.sigma, size=self.current_geometry.shape)
        old_accelerations = np.copy(self.current_accelerations)

        # update geometry, forces, and velocities
        self.current_geometry += self.b[:, np.newaxis] * self.dt * (self.current_velocities + \
            0.5 * self.dt *  self.current_accelerations + \
            0.5 * gaussian_random_numbers / self.masses[:, np.newaxis])
        self.update_energy_and_accelerations()
        self.current_velocities = self.a[:, np.newaxis] * self.current_velocities + \
            0.5 * self.dt * (self.a[:, np.newaxis] * old_accelerations + self.current_accelerations) \
            + (self.b / self.masses)[:, np.newaxis] * gaussian_random_numbers
        self.center_positions()

################################################################################################

class Langevin_Thermostat(Thermostat):
    """
    Implements the Langevin theromstat which models the system via Langevin dynamics. This is
    a stochastic thermostat, so is not time-reversible, but correctly samples the canonical 
    ensemble and also provides completely valid dynamics, so that time correlation functions 
    can be calculated from the trajectory. A single friction parameter is provided.

    This is adapted from eqs. 29 & 30 of:
    Farago, Langevin thermostat for robust configurational and kinetic sampling
    This is the "GJF-2GJ" algorithm.
    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature: float,
                alpha=25.0):
        super().__init__(geometry, masses, potential_function, dt, temperature)
        self.alpha = alpha

        self.a = (1 - self.alpha * self.dt / (2 * self.masses)) / (1 + self.alpha * self.dt / (2 * self.masses))
        self.sqrtb = np.sqrt(1 / (1 + self.alpha * self.dt / (2 * self.masses)))
        self.sigma = np.sqrt(2 * self.alpha * self.temperature * self.dt)
        self.gaussian_noise = np.random.default_rng().normal(0.0, self.sigma, size=self.current_geometry.shape)
    
    def integrate(self):
        # get the random numbers and copies of previous random numbers from previous step
        old_gaussian_noise = np.copy(self.gaussian_noise)
        self.gaussian_noise = np.random.default_rng().normal(0.0, self.sigma, size=self.current_geometry.shape)

        # update geometry, forces, and velocities
        self.current_velocities = self.a[:, np.newaxis] * self.current_velocities + \
            self.sqrtb[:, np.newaxis] * self.dt * self.current_accelerations + \
            (self.sqrtb / (2 * self.masses))[:, np.newaxis] * (self.gaussian_noise + old_gaussian_noise)
        
        self.current_geometry += self.sqrtb[:, np.newaxis] * self.current_velocities * self.dt

        self.update_energy_and_accelerations()
        self.center_positions()
