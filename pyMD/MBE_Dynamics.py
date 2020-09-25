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
    mbe_ff = MBE_Potential(4, fragments, ttm21f)

    geometry = np.vstack(fragments.fragments)
    
    atom_masses = np.array(list(map(get_mass_of_element, fragments.atom_labels)))
    

    equilibration_output = [
        ("step",                "w6_dynamics_equil_full.md"),
        ("time",                "w6_dynamics_equil_full.md"),
        ("potential_energy",    "w6_dynamics_equil_full.md"),
        ("kinetic_energy",      "w6_dynamics_equil_full.md"),
        ("total_energy",        "w6_dynamics_equil_full.md"),
        ("temperature",         "w6_dynamics_equil_full.md"),
        ("geometry",            "w6_dynamics_geometry_equil_full.xyz"),
        ("velocity",            "w6_dynamics_velocity_equil_full.xyz")
    ]

    production_output = [
        ("step",                "w6_dynamics_production_full.md"),
        ("time",                "w6_dynamics_production_full.md"),
        ("potential_energy",    "w6_dynamics_production_full.md"),
        ("kinetic_energy",      "w6_dynamics_production_full.md"),
        ("total_energy",        "w6_dynamics_production_full.md"),
        ("temperature",         "w6_dynamics_production_full.md"),
        ("geometry",            "w6_dynamics_geometry_production_full.xyz"),
        ("velocity",            "w6_dynamics_velocity_production_full.xyz"),
        ("force",               "w6_dynamics_forces_production_full.xyz")
    ]

    # time step
    ts=20.7

    #mbe_ff.evaluate_on_geometry,
    ####### INTEGRATORS FOR EQUILIBRATION AND PRODUCTION #######
    thermostatted_equilibration = Langevin_Thermostat(geometry, 
                                     atom_masses, 
                                     ttm21f.evaluate,
                                     dt=ts,
                                     temperature=300,
                                     alpha=25.0)

    nve_production = Velocity_Verlet(geometry, 
                                     atom_masses, 
                                     ttm21f.evaluate,
                                     dt=ts,
                                     initial_velocities=initial_velocities)

    ####### LOGGERS FOR EQUILIBRATION AND PRODUCTION #######
    log_equil      = Logger(equilibration_output, 5, fragments.atom_labels)
    log_production = Logger(production_output, 5, fragments.atom_labels)
    
    ####### DYNAMICS MANAGERS FOR EQUILIBRATION AND PRODUCTION #######
    dynamics_equil =      Dynamics(thermostatted_equilibration,
                                   log_equil,
                                   max_iterations=100000) # .05 ns equilibration in NVT

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