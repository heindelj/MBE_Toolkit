import numpy as np
from types import FunctionType

class Integrator:
    """
    Abstract base class from which all real integration methods will inherit.
    This includes any thermostats, which are themselves integrators of the equations-
    of-motion.
    Args:
        geometry            : initial geometry to begin dynamics from
        masses              : masses of all atoms given in geometry
        potential           : potential function (likely derived from a potential object)
        dt                  : time step for the integration (in atomic units)
        initital_velocities : initial velocities from which to begin

    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature=None,
                initial_velocities=None):

        self.current_geometry = geometry * 1.88973 # convert to bohr UNITS
        self.masses = masses * 1822.888486 # convert to atomic unit of masses UNITS
        self.potential_function = potential_function
        self.extras_from_potential = {} # extras dictionary which may or may not passed from potential
        self.dt = dt

        self.current_velocities = np.zeros((len(self.masses), 3))

        if initial_velocities is not None:
            self.current_velocities = initial_velocities

        # initialize velocities from maxwell-boltzmann distribution
        if temperature is not None:
            self.temperature = temperature / 315775.0248 # convert to atomic unit of temperature UNITS
            self.current_velocities = self.sample_maxwell_boltzmann_velocities()
            self.rescale_velocities()

        # initialize energy and forces based on initial position
        # there is an argument to be made that the forces should be initialized to zero
        # but it shouldn't make much of a difference either way
        self.update_energy_and_accelerations()
    
    def integrate(self):
        raise NotImplementedError("You tried to integrate with the base class integrator. This is an abstract base class. Please choose a child of this class instead.")
    
    def update_energy_and_accelerations(self):
        """Makes the call to potential function"""
        potential_data = self.potential_function(self.current_geometry / 1.88973) #UNITS
        if len(potential_data) == 2:
            self.energy, self.current_accelerations = potential_data
        elif len(potential_data) == 3:
            self.energy, self.current_accelerations, self.extras_from_potential = potential_data
        self.current_accelerations /= self.masses[:, np.newaxis]

    def sample_maxwell_boltzmann_velocities(self, indices=None):
        """
        Samples a Maxwell-Boltzmann distribution for velocities.
        Args:
            indices: a list of indices to sample instead of all atoms.
        Returns:
            velocities: numpy array of velocities sampled from M-B and rescaled to match target temperature
        """
        if self.temperature:
            if indices:
                sigmas = np.sqrt(self.temperature * (1.0 / self.masses[indices]))
                velocities = np.zeros((len(indices), 3))
            else:
                sigmas = np.sqrt(self.temperature * (1.0 / self.masses))
                velocities = np.zeros((len(self.masses), 3))
            for i, sigma in enumerate(sigmas):
                velocities[i,:] = np.random.default_rng().normal(0.0, sigma, size=3)
            return velocities
        else:
            raise ValueError("There is no temperature set. We cannot sample the velocities without a temperature. Check that you have provided a temperature to the default integrator.")

    def rescale_velocities(self):
        self.current_velocities *= np.sqrt(self.temperature / self.get_temperature())

    def get_kinetic_energy(self):
        kinetic_energy = 0.5 * self.masses * np.einsum('ij,ij->i', self.current_velocities, self.current_velocities)
        return np.sum(kinetic_energy)

    def get_momentum(self):
        return self.current_velocities * self.masses[:, np.newaxis]
    
    def get_temperature(self):
        # TODO: Make DOFs a parameter which is determined at construction based on
        # an option which is given for removing COM motion and rotational motion
        return 2 * self.get_kinetic_energy() / (3 * len(self.masses) )# - 3)

    def remove_com_motion(self):
        mean = self.current_velocities.mean(axis=0) # create a view instead of copy
        self.current_velocities -= mean
    
    def center_positions(self):
        mean = self.current_geometry.mean(axis=0) # create a view instead of copy
        self.current_geometry -= mean


###############################################################################################

class Verlet(Integrator):
    """
    Traditional Verlet integrator.
    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature=None,
                initial_velocities=None):
        super().__init__(geometry, masses, potential_function, dt, temperature=temperature, initial_velocities=initial_velocities)
        self.old_geometry = np.copy(self.current_geometry)
    
    def integrate(self):
        """
        Implements verlet integration with midpoint average for the velocity
        """
        # propagate geometry
        new_geometry = 2 * self.current_geometry - self.old_geometry + self.current_accelerations * self.dt**2

        # get velocity from midpoint of x(t+1) and x(t-1)
        self.current_velocities = (new_geometry - self.old_geometry) / (2 * self.dt)

        # update stored geometries
        self.old_geometry     = np.copy(self.current_geometry)
        self.current_geometry = np.copy(new_geometry)

        self.update_energy_and_accelerations()
        self.center_positions()

###############################################################################################

class Velocity_Verlet(Integrator):
    """
    Velocity-Verlet integrator.
    """
    def __init__(self,
                geometry: np.ndarray,
                masses: np.ndarray,
                potential_function: FunctionType,
                dt: float,
                temperature=None,
                initial_velocities=None):
        super().__init__(geometry, masses, potential_function, dt, temperature=temperature, initial_velocities=initial_velocities)

    def integrate(self):
        """
        Implements the velocity verlet algorithm.
        """
        # half-step of velocities
        self.current_velocities += 0.5 * self.current_accelerations * self.dt

        # update geometry
        self.current_geometry += self.current_velocities * self.dt
        
        self.update_energy_and_accelerations()

        # second half-step of velocities
        self.current_velocities += 0.5 * self.current_accelerations * self.dt
        self.center_positions()