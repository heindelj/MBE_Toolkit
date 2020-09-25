
import sys
import numpy as np

__all__ = ['Potential', 'TTM']

class Potential:
    """Abstract base class for potential energy surfaces. A single function, evaluate(),
    must be implemented which returns the energy and gradients.

    Each child class of this needs to implement any methods to:
    1) returns the gradients with the same ordering as it received them
    2) help with evaluation of the potential
    """
    def __init__(self):
        pass
    
    def evaluate(self, coords):
        raise NotImplementedError

class TTM(Potential):
    import os
    lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../bin")
    os.environ['LD_LIBRARY_PATH'] = lib_path
    try:
        import ttm
    except ImportError:
        print("Did not find ttm module. Make sure you have compiled it and the TTM library can be linked against.")
        sys.exit(1)
    def __init__(self, model=21):
        """Evaluates the energy and gradients of the TTM family of potentials.

        Args:
            model (int, optional): The TTM model which will be used. Options are 2, 21, and 3. Defaults to 21.
        """
        self.model = model
        possible_models = [2, 21, 3]
        if self.model not in possible_models:
            print("The possible TTM versions are 2, 21, or 3. Please choose one of these.")
            sys.exit(1)
    
    def evaluate(self, coords):
        """Takes xyz coordinates of water molecules in O H H, O H H order and re-orders to OOHHHH order
        then transposes to fortran column-ordered matrix and calls the TTM potential from an f2py module.


        Args:
            coords (ndarray3d): xyz coordinates of a system which can be evaluated by this potential.
        Returns:
            energy (float): energy of the system in hartree
            forces (ndarray3d): forces of the system in hartree / bohr
        """
        #Sadly, we need to re-order the geometry to TTM format which is all oxygens first.
        coords = self.ttm_ordering(coords)
        gradients, energy = self.ttm.ttm_from_f2py(self.model, np.asarray(coords).T, int(len(coords) / 3))
        return energy / 627.5, (-self.normal_water_ordering(gradients.T) / 627.5) / 1.88973
    
    @staticmethod
    def ttm_ordering(coords):
        """Sorts an array of coordinates in OHHOHH format to OOHHHH format.

        Args:
            coords (ndarray3d): numpy array of coordinates

        Returns:
            ndarray3d: numpy array of coordinate sorted according to the order TTM wants.
        """
        atom_order = []
        for i in range(0, coords.shape[0], 3):
            atom_order.append(i)
        for i in range(0, coords.shape[0], 3):
            atom_order.append(i+1)
            atom_order.append(i+2)
        return coords[atom_order,:]
    
    @staticmethod
    def normal_water_ordering(coords):
        """Sorts an array of coordinates in OOHHHH format to OHHOHH format.

        Args:
            coords (ndarray3d): numpy array of coordinates

        Returns:
            ndarray3d: numpy array of coordinate sorted in the normal way for water.
        """
        atom_order = []
        Nw = int(coords.shape[0] / 3)
        for i in range(0, Nw, 1):
            atom_order.append(i)
            atom_order.append(Nw+2*i)
            atom_order.append(Nw+2*i+1)
        return coords[atom_order,:]