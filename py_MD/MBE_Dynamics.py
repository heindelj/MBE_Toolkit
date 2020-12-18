from .Dynamics import Dynamics
from .MBE_Potential import MBE_Potential
from .Fragments import Fragments
from .Potential import *
from .Integrator import *
from .Thermostat import *
from .Masses import get_mass_of_element
from .Logger import Logger

import numpy as np
import sys, os, time


if __name__ == '__main__':
    try:
        ifile, mbe_order, input_temperature = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    except:
        print("[xyz file] [mbe order] [temperature]")
        sys.exit(1)
    
    fragments = Fragments(ifile)
    ttm21f = TTM("/home/heindelj/Research/Sotiris/MBE_Dynamics/MBE_Dynamics_Home_Code/pyMD/bin")
    mbe_ff = MBE_Potential(mbe_order, fragments, ttm21f, nproc=8, return_extras=True)

    geometry = np.vstack(fragments.fragments)
    
    atom_masses = np.array(list(map(get_mass_of_element, fragments.atom_labels)))
    

    equilibration_output = [
        ("step",                              "w10_dynamics_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("time",                              "w10_dynamics_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("potential_energy",                  "w10_dynamics_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("kinetic_energy",                    "w10_dynamics_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("total_energy",                      "w10_dynamics_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("temperature",                       "w10_dynamics_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("geometry",                 "w10_dynamics_geometry_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.xyz"),
        ("velocity",                 "w10_dynamics_velocity_equil_"+str(input_temperature)+"K_"+str(mbe_order)+"body.xyz"),
        ("extras", None)
    ]

    production_output = [
        ("step",                         "w10_dynamics_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("time",                         "w10_dynamics_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("potential_energy",             "w10_dynamics_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("kinetic_energy",               "w10_dynamics_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("total_energy",                 "w10_dynamics_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("temperature",                  "w10_dynamics_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.md"),
        ("geometry",            "w10_dynamics_geometry_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.xyz"),
        ("velocity",            "w10_dynamics_velocity_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.xyz"),
        ("force",                 "w10_dynamics_forces_production_"+str(input_temperature)+"K_"+str(mbe_order)+"body.xyz")
    ]

    # time step
    ts=20.7

    #mbe_ff.evaluate_on_geometry
    ####### INTEGRATORS FOR EQUILIBRATION AND PRODUCTION #######
    thermostatted_equilibration = Langevin_Thermostat(geometry, 
                                     atom_masses, 
                                     mbe_ff.evaluate_on_geometry_parallel,
                                     dt=ts,
                                     temperature=input_temperature,
                                     alpha=25.0)

    nve_production = Velocity_Verlet(geometry,
                                     atom_masses,
                                     mbe_ff.evaluate_on_geometry_parallel,
                                     temperature=input_temperature,
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
    dynamics_production.integrator.rescale_velocities()
    dynamics_production.integrator.current_geometry = dynamics_equil.integrator.current_geometry
    dynamics_production.integrator.current_accelerations = dynamics_equil.integrator.current_accelerations
    dynamics_production.propagate()
    end = time.time()
    print("Equilibration Time: ", end - start, " seconds")