from numpy.lib.function_base import gradient
from Fragments import Fragments
from Potential import *
from MBE_Potential import MBE_Potential

class Composite_Potential:
    """
    A composition of multiple Potential objects which are used to construct MBE_Potentials
    that are used to calculate all orders of the MBE.
    """
    def __init__(self, orders_and_potentials: dict, fragments: Fragments, full_background_potential=True):
        """
        Takes a dictionary of integers specifying the maximum order of the MBE and corresponding potential which
        will be used for this method.
        """
        self.orders_and_potentials = orders_and_potentials
        self.fragments = fragments
        self.full_background_potential = full_background_potential
        self.mbe_potentials = []
        self.create_member_potentials()

    def create_member_potentials(self, parallel=False):
        # get the maximum order which will subtract out 
        # the residual from full minus MBE up to the max order
        max_order = max(self.orders_and_potentials.keys())
        max_order_potential = self.orders_and_potentials[max_order]
        
        # get the minimum order which should return 
        # the sum up to minimum order
        min_order = min(self.orders_and_potentials.keys())
        self.mbe_potentials.append(MBE_Potential(min_order, self.fragments, self.orders_and_potentials[min_order]))
        for mbe_order, mbe_potential in self.orders_and_potentials.items():
            if mbe_order != max_order and mbe_order != min_order:
                self.mbe_potentials.append(MBE_Potential(mbe_order, self.fragments, mbe_potential, return_order_n=mbe_order))
        # First potential does MBE up to max_order - 1. Second potential does the full calculation and then subtracts out the difference
        self.full_system_potential = (MBE_Potential(max_order-1, self.fragments, max_order_potential), max_order_potential)

    def get_energy_and_gradients(self, coords, parallel_MBE=False):
        """
        Calls all of the potentials contained in self.mbe_potentials and self.full_system_potential
        to get the total energy for this composite potential.
        """
        nbody_energies = np.zeros(len(self.mbe_potentials)+1)
        total_gradients = np.zeros_like(coords)
        for (i, potential) in enumerate(self.mbe_potentials):
            if parallel_MBE:
                energy, gradients = potential.evaluate_on_geometry_parallel(coords)
            else:
                energy, gradients = potential.evaluate_on_geometry(coords)
            nbody_energies[i] = energy
            total_gradients += gradients
        
        residual_energy_mbe, residual_gradients_mbe = self.full_system_potential[0].evaluate_on_geometry(coords)
        residual_energy_full, residual_gradients_full = self.full_system_potential[1].evaluate(coords)
        nbody_energies[-1] = residual_energy_full - residual_energy_mbe
        total_gradients   += residual_gradients_full - residual_gradients_mbe
        return np.sum(nbody_energies), total_gradients
    
if __name__ == '__main__':
    import sys
    try:
        ifile = sys.argv[1]
    except:
        print("Didn't get an xyz file.")
        sys.exit(1)
    
    import numpy as np
    import time

    fragments = Fragments(ifile)
    ttm21f = TTM("/home/heindelj/dev/python_development/MBE_Toolkit/bin/")
    mbpol = MBPol("/home/heindelj/dev/python_development/MBE_Toolkit/bin/")

    orders_and_potentials = {3: ttm21f, 2: mbpol}

    comp_potential = Composite_Potential(orders_and_potentials, fragments)
    start = time.time()
    energy, gradients = comp_potential.get_energy_and_gradients(np.vstack(fragments.fragments), parallel_MBE=True)
    print(gradients)
    print("Total Energy Composite MBE: ", "{:.6f}".format(energy * 627.5), "kcal/mol")
    print(time.time() - start, " seconds")


    #ttm_mbe = MBE_Potential(3, fragments, ttm21f)
    #mbpol_mbe = MBE_Potential(2, fragments, mbpol)
