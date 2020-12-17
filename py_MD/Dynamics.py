from .Logger import Logger
from .Integrator import Integrator

from typing import Callable
import numpy as np
import sys

class Dynamics:
    """Manages the dynamics simulation by making calls to the integrator which does the
    work of the propagation and makes calls to the potential. Also implements all of the
    logging methods and calls the logger object which manages which files data is written
    to and the strides of when to write, etc.

    TODO: Long-Term: Dynamics should take an ensemble class which manages the Thermostats and
    Barostats. The Integrator should probably get passed into this as well?
    """
    def __init__(self,
                 integrator: Integrator,
                 log: Logger,
                 max_iterations: int):
        
        self.integrator = integrator
        self.log = log

        self.max_iterations = max_iterations
        self.current_step = 1

        # logging dictionary
        # all of these need to be update to grab data from the integrator.
        self.logging_vtable = {
            "step":                self.log_step,
            "time":                self.log_time,
            "momentum":            self.log_momentum,
            "potential_energy":    self.log_potential_energy,
            "kinetic_energy":      self.log_kinetic_energy,
            "mean_kinetic_energy": self.log_average_kinetic_energy,
            "total_energy":        self.log_total_energy,
            "temperature":         self.log_temperature_from_kinetic_energy,
            "geometry":            self.log_geometry,
            "velocity":            self.log_velocity,
            "force":               self.log_forces,
            "extras":              self.log_extras
        }
    

    def propagate(self):
        """Propagates current geometry forward in time using self.potential_function
        via velocity-verlet integration.
        """
        for self.current_step in range(self.max_iterations):
            self.integrator.integrate()
            self.log_data()

    def log_data(self):
        """Calls the getter function stored in self.logging_vtable for a pre-defined set of
        observables which the user can ask for via self.log. If the observable is not 
        in the dictionary, the request is ignored.

        Passing the key prevents possible typos between the vtable and getter function.
        """
        if (self.get_step() % self.log.logging_stride) == 0:
            for key in self.log.keys:
                if key in self.logging_vtable:
                    self.logging_vtable[key](key)
                else:
                    print(f"Received invalid observable. Please check that there is a getter for the observable {key}")
                    sys.exit(1)
        self.log.log()

    def get_step(self):
        return self.current_step

    def log_step(self, key):
        self.log.data_log[key].append(self.current_step)

    def log_time(self, key):
        self.log.data_log[key].append(self.current_step * self.integrator.dt)

    def log_geometry(self, key):
        self.log.data_log[key].append(self.integrator.current_geometry / 1.88973) # bohr to angstrom UNITS
    
    def log_velocity(self, key):
        self.log.data_log[key].append(self.integrator.current_velocities / 1.88973) # bohr to angstrom UNITS

    def log_forces(self, key):
        self.log.data_log[key].append(self.integrator.current_accelerations * self.integrator.masses[:, np.newaxis])

    def log_extras(self, key=None):
        for key, value in self.integrator.extras_from_potential.items():
            self.log.data_log[key].append(value)
            if key not in self.log.file_log:
                self.log.file_log[key] = key + ".dat"
                # just open to delete previous file
                with open(self.log.file_log[key], 'w'):
                    pass
                self.log.special_logging_keys.append(key)

    def log_temperature_from_kinetic_energy(self, key):
        """Returns the temperature according to the average kinetic energy in kelvin.
        <E>=3/2kT; T=2/3<E> (k=1); <E>=1/2m<v^2> where <v^2> is the rms velocity
        """
        self.log.data_log[key].append(self.integrator.get_temperature()  * 315775.0248) #UNITS
    
    def log_momentum(self, key):
        self.log.data_log[key].append(self.integrator.get_momentum())

    def log_potential_energy(self, key):
        self.log.data_log[key].append(self.integrator.energy)

    def log_kinetic_energy(self, key):
        self.log.data_log[key].append(self.integrator.get_kinetic_energy())
    
    def log_average_kinetic_energy(self, key):
        momentum = self.integrator.get_momentum()
        kinetic_energy = np.einsum('ij,ij->i', momentum, momentum) / (2 * self.integrator.masses)
        self.log.data_log[key].append(np.mean(kinetic_energy))

    def log_total_energy(self, key):
        self.log.data_log[key].append(self.integrator.energy + self.integrator.get_kinetic_energy())