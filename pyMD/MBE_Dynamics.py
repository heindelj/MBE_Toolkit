from Dynamics import Dynamics
from MBE_Potential import MBE_Potential
from Fragments import Fragments
from Potential import *
from Integrator import *
from Thermostat import *
from Masses import get_mass_of_element
from Logger import Logger

import numpy as np
import sys, time


if __name__ == '__main__':
    try:
        ifile = sys.argv[1]
    except:
        print("Didn't get an xyz file.")
        sys.exit(1)
    
    fragments = Fragments(ifile)
    ttm21f = TTM(21)
    mbe_ff = MBE_Potential(2, fragments, ttm21f)

    geometry = np.vstack(fragments.fragments)
    
    atom_masses = np.array(list(map(get_mass_of_element, fragments.atom_labels)))
    

    equilibration_output = [
        ("step",                              "w6_dynamics_equil_2body.md"),
        ("time",                              "w6_dynamics_equil_2body.md"),
        ("potential_energy",                  "w6_dynamics_equil_2body.md"),
        ("kinetic_energy",                    "w6_dynamics_equil_2body.md"),
        ("total_energy",                      "w6_dynamics_equil_2body.md"),
        ("temperature",                       "w6_dynamics_equil_2body.md"),
        ("geometry",                 "w6_dynamics_geometry_equil_2body.xyz"),
        ("velocity",                 "w6_dynamics_velocity_equil_2body.xyz")
    ]

    production_output = [
        ("step",                         "w6_dynamics_production_2body.md"),
        ("time",                         "w6_dynamics_production_2body.md"),
        ("potential_energy",             "w6_dynamics_production_2body.md"),
        ("kinetic_energy",               "w6_dynamics_production_2body.md"),
        ("total_energy",                 "w6_dynamics_production_2body.md"),
        ("temperature",                  "w6_dynamics_production_2body.md"),
        ("geometry",            "w6_dynamics_geometry_production_2body.xyz"),
        ("velocity",            "w6_dynamics_velocity_production_2body.xyz"),
        ("force",                 "w6_dynamics_forces_production_2body.xyz")
    ]

    # time step
    ts=20.7

    ####### INTEGRATORS FOR EQUILIBRATION AND PRODUCTION #######
    thermostatted_equilibration = Langevin_Thermostat(geometry, 
                                     atom_masses, 
                                     mbe_ff.evaluate_on_geometry,
                                     dt=ts,
                                     temperature=300,
                                     alpha=25.0)

    nve_production = Velocity_Verlet(geometry, 
                                     atom_masses, 
                                     mbe_ff.evaluate_on_geometry,
                                     dt=ts)

    ####### LOGGERS FOR EQUILIBRATION AND PRODUCTION #######
    log_equil      = Logger(equilibration_output, 5, fragments.atom_labels)
    log_production = Logger(production_output, 5, fragments.atom_labels)
    
    ####### DYNAMICS MANAGERS FOR EQUILIBRATION AND PRODUCTION #######
    dynamics_equil =      Dynamics(thermostatted_equilibration,
                                   log_equil,
                                   max_iterations=200000) # .1 ns equilibration in NVT

    dynamics_production = Dynamics(nve_production,
                                   log_production,
                                   max_iterations=2000000) # 1 ns production run in NVE
    start = time.time()
    dynamics_equil.propagate()
    end = time.time()
    print("Equilibration Time: ", end - start, " seconds")
    start = time.time()
    # TODO Implement class method to initialize one Integrator from another Integrators state
    dynamics_production.integrator.current_velocities = dynamics_equil.integrator.current_velocities
    dynamics_production.integrator.current_geometry = dynamics_equil.integrator.current_geometry
    dynamics_production.integrator.current_accelerations = dynamics_equil.integrator.current_accelerations
    dynamics_production.propagate()
    end = time.time()
    print("Equilibration Time: ", end - start, " seconds")