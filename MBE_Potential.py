from Fragments import Fragments
from Dynamics import Dynamics
from Potential import *
import numpy as np
from math import comb
import sys, os, time
import multiprocessing as mp

class MBE_Potential:
    """
    Implements an MBE potential which calls out to a Potential object and
    parses the output energy and forces
    """
    def __init__(self, highest_order: int, fragments: Fragments, potential: Potential):
        self.highest_order = highest_order
        self.fragments = fragments
        self.potential = potential
        nproc = 8
        self._pool = mp.Pool(nproc)
    
    def evaluate_on_fragments(self):
        """
        Uses the Potential object passed in to calculate the forces and energies for every fragment 
        in self.fragments.

        This operates directly on the fragments brought in with self.fragments
        """
        # this will hang on to the n-body forces and energies summed for all
        # n-mers. These will then be combined with proper combinatorial weights
        # to get the total forces and energies.
        forces_sum = []
        energy_sum = []

        for order in range(self.highest_order):
            fragment_combinations = self.fragments.make_nmers(order + 1)
            atom_indices = self.fragments.get_indices_for_fragment_combination(order + 1)

            # each force sum is the size of forces for total system
            forces_sum.append(np.zeros( (len(self.fragments.atom_labels), 3), dtype=np.float64))
            energy_sum.append(0.0)

            # now get components of the MBE by summing energies and forces for each order
            for i_frag, fragment in enumerate(fragment_combinations):
                energy, forces = self.potential.evaluate(fragment)

                # add forces in the appropriate rows of force sum
                energy_sum[order] += energy
                assert(len(forces) == len(fragment))
                for i_force, force in enumerate(forces):
                    forces_sum[order][atom_indices[i_frag][i_force]] += force

        ### now that we have the sums for each part, we must weight them by the 
        ### combinatorial number of times they show up and accumulate the totals.
        total_forces = np.zeros( (len(self.fragments.atom_labels), 3), dtype=np.float64)
        total_energy = 0.0

        nbody_energies = [0 for x in range(self.highest_order)]
        nbody_forces = [np.zeros_like(total_forces) for x in range(self.highest_order)]
        
        # set the appropriate 1-body terms
        nbody_forces[0] = forces_sum[0]
        nbody_energies[0] = energy_sum[0]
        
        # multiply each many-body sum by the appropriate combinatorial factor
        N = len(self.fragments.fragments)
        for iMBE in range(1, self.highest_order):
            for i in range(iMBE+1):
                nbody_energies[iMBE] += (-1)**i * comb(N-(iMBE+1)+i,i) * energy_sum[iMBE-i]
                nbody_forces[iMBE]   += (-1)**i * comb(N-(iMBE+1)+i,i) * forces_sum[iMBE-i]

        # accumulate many-body energies and forces
        total_energy = np.sum(nbody_energies)
        total_forces = np.sum(nbody_forces, axis=0)

        #print(total_energy * 627.5)
        #print(total_forces * 627.5 / 1.88973)

        return total_energy, total_forces

    def evaluate_on_fragments_parallel(self):
        """
        Uses the Potential object passed in to calculate the forces and energies for every fragment 
        in self.fragments.

        This operates directly on the fragments brought in with self.fragments
        """
        # this will hang on to the n-body forces and energies summed for all
        # n-mers. These will then be combined with proper combinatorial weights
        # to get the total forces and energies.
        all_fragments = []
        all_indices_into_fragments = []
        for order in range(self.highest_order):
            all_fragments += self.fragments.make_nmers(order + 1)
            all_indices_into_fragments += self.fragments.get_indices_for_fragment_combination(order + 1)

        # evaluate the potential on all fragments
        energies, forces = map(list, zip(*self._pool.map(self.potential.evaluate, all_fragments)))

        # each force sum is the size of forces for total system
        forces_sum = np.zeros( (self.highest_order, len(self.fragments.atom_labels), 3), dtype=np.float64)
        energy_sum = np.zeros( (self.highest_order,) )

        # loop through the appropriate parts of the energies and forces arrays and sum them for each n-body term
        N = len(self.fragments.fragments)
        comb_sum = 0
        global_fragment_index = 0
        for order in range(self.highest_order):
            mask = np.full(len(energies), False)
            mask[np.arange(comb_sum, comb(N, order+1) + comb_sum)] = True
            energy_sum[order] = np.sum(energies, where=mask)

            for nbody_forces in forces[comb_sum:(comb(N, order+1) + comb_sum)]:
                indices = all_indices_into_fragments[global_fragment_index]
                global_fragment_index += 1
                for i_force, force in enumerate(nbody_forces):
                    forces_sum[order][indices[i_force]] += force
            
            comb_sum += comb(N, order+1)

        ### now that we have the sums for each part, we must weight them by the 
        ### combinatorial number of times they show up and accumulate the totals.
        total_forces = np.zeros( (len(self.fragments.atom_labels), 3), dtype=np.float64)
        total_energy = 0.0

        nbody_energies = [0 for x in range(self.highest_order)]
        nbody_forces = [np.zeros_like(total_forces) for x in range(self.highest_order)]
        
        # set the appropriate 1-body terms
        nbody_forces[0] = forces_sum[0]
        nbody_energies[0] = energy_sum[0]
        
        # multiply each many-body sum by the appropriate combinatorial factor
        N = len(self.fragments.fragments)
        for iMBE in range(1, self.highest_order):
            for i in range(iMBE+1):
                nbody_energies[iMBE] += (-1)**i * comb(N-(iMBE+1)+i,i) * energy_sum[iMBE-i]
                nbody_forces[iMBE]   += (-1)**i * comb(N-(iMBE+1)+i,i) * forces_sum[iMBE-i]

        # accumulate many-body energies and forces
        total_energy = np.sum(nbody_energies)
        total_forces = np.sum(nbody_forces, axis=0)

        #print(total_energy * 627.5)
        #print(total_forces * 627.5 / 1.88973)

        return total_energy, total_forces

    def evaluate_on_geometry(self, geometry):
        """This is a thin wrapper around evaluate_on_fragments() which allows
        raw coordinates to be passed in, and then fragments those coordinates
        according to the shape of the self.fragments.fragments.

        Args:
            geometry (ndarray): Nx3 array of cartesian coordinates
        """
        self.fragments.fragment_geometry(geometry)
        energy, forces = self.evaluate_on_fragments()
        return energy, forces
    
    def evaluate_on_geometry_parallel(self, geometry):
        """This is a thin wrapper around evaluate_on_fragments() which allows
        raw coordinates to be passed in, and then fragments those coordinates
        according to the shape of the self.fragments.fragments.

        Args:
            geometry (ndarray): Nx3 array of cartesian coordinates
        """
        self.fragments.fragment_geometry(geometry)
        energy, forces = self.evaluate_on_fragments_parallel()
        return energy, forces

if __name__ == '__main__':
    try:
        ifile = sys.argv[1]
    except:
        print("Didn't get an xyz file.")
        sys.exit(1)
    
    fragments = Fragments(ifile)
    ttm21f = TTM(21)
    mbe_ff = MBE_Potential(4, fragments, ttm21f)
    
    start = time.time()
    mbe_ff.evaluate_on_fragments_parallel()
    print(time.time() - start)
