import numpy as np
from Logger import Logger
import Integrator

from typing import Callable
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
            "step":                self.get_step,
            "time":                self.get_time,
            "momentum":            self.get_momentum,
            "potential_energy":    self.get_potential_energy,
            "kinetic_energy":      self.get_kinetic_energy,
            "mean_kinetic_energy": self.get_average_kinetic_energy,
            "total_energy":        self.get_total_energy,
            "temperature":         self.get_temperature_from_kinetic_energy,
            "geometry":            self.get_geometry,
            "velocity":            self.get_velocity,
            "force":               self.get_forces
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

        Refactor: This needs to be here because of the call to the logging_vtable.
        Maybe there is a way to allow the logging class to call the functions in the vtable?
        """
        if (self.get_step() % self.log.logging_stride) == 0:
            for key in self.log.keys:
                if key in self.logging_vtable:
                    self.log.data_log[key].append(self.logging_vtable[key]())
                else:
                    print(f"Received invalid observable. Please check that there is a getter for the observable {key}")
                    sys.exit(1)
        self.log.log()

    def get_step(self):
        return (self.current_step)

    def get_time(self):
        return (self.current_step) * self.integrator.dt

    def get_geometry(self):
        return self.integrator.current_geometry / 1.88973 # bohr to angstrom UNITS
    
    def get_velocity(self):
        return self.integrator.current_velocities / 1.88973 # bohr to angstrom UNITS

    def get_forces(self):
        return self.integrator.current_accelerations * self.integrator.masses[:, np.newaxis]

    def get_temperature_from_kinetic_energy(self):
        """Returns the temperature according to the average kinetic energy in kelvin.
        <E>=3/2kT; T=2/3<E> (k=1); <E>=1/2m<v^2> where <v^2> is the rms velocity
        """
        return self.integrator.get_temperature()  * 315775.0248 #UNITS
    
    def get_momentum(self):
        return self.integrator.current_velocities * self.integrator.masses[:, np.newaxis]

    def get_potential_energy(self):
        return self.integrator.current_energy

    def get_kinetic_energy(self):
        return self.integrator.get_kinetic_energy()
    
    def get_average_kinetic_energy(self):
        momentum = self.get_momentum()
        kinetic_energy = np.einsum('ij,ij->i', momentum, momentum) / (2 * self.integrator.masses)
        return np.mean(kinetic_energy)

    def get_total_energy(self):
        return self.integrator.current_energy + self.get_kinetic_energy()