import sys, os, subprocess, importlib
import numpy as np
import ctypes

__all__ = ['Potential', 'TTM', 'MBPol']

class Potential:
    """Abstract base class for potential energy surfaces. A single function, evaluate(),
    must be implemented which returns the energy and gradients.

    Each child class needs to implement any methods to:
    1) returns the gradients with the same atom ordering as it received them
    2) help with evaluation of the potential by calling a library or python function
    """
    def __init__(self, path_to_library=None, name_of_function=None, name_of_module=None, name_of_library=None):
        """
        path_to_library  (str): absolute path to library or file which contains potential function.
        name_of_function (str): name of the function to be called from the imported module.
        name_of_module   (str): name of a module containing the function to be called.
        name_of_library  (str): name of a shared library containing the function to be called.
        """
        self.path_to_library = os.path.normpath(os.path.join(os.getcwd(), path_to_library))
        self.name_of_function = name_of_function
        self.name_of_module = name_of_module
        self.name_of_library = name_of_library
        self.work_dir = os.getcwd()
        
        # these are the actual handles to the library, module, and/or function
        self.library = None
        self.potential_function = None

        if self.name_of_module and self.name_of_library:
            print("Please provide either a module name or a library name, but not both.")
            print("name_of_module is for a python module we will call.")
            print("name_of_library is for a shared a library we will call into.")
            sys.exit(1)

        # remove the .py suffix from name_of_module if it happens to be there
        if self.name_of_module and ".py" in self.name_of_module:
            self.name_of_module = self.name_of_module.split(".")[0]
    
    def evaluate(self, coords):
        raise NotImplementedError

    def initialize_potential(self):
        """
        Initializes a potential which is accessed via an absolute path, self.path_to_library, and a function name.
        Either a module name or the name of a shared library may be provided for python functions and C/Fortran functions, respectively.

        This function is crucial if you want to use a multiprocessing pool, as it needs to be
        used as the initializer for the pool.
        """
        if self.potential_function is None:
            os.chdir(self.path_to_library)
            sys.path.insert(0, os.getcwd())
            # this branch is for loading a function from a python module
            if self.name_of_module:
                try:
                    module = importlib.import_module(self.name_of_module)
                    self.potential_function = getattr(module, self.name_of_function)
                    os.chdir(self.work_dir)
                except ImportError:
                    print("Did not find potential module. Make sure you have compiled it and the library can be linked against, including things like libgfortran and libgcc.")
                    print("If the module is a plain python function, then make sure you are passing the correct absolute path to the file.")
                    sys.exit(1)
            elif self.name_of_library:
                try:
                    self.library = ctypes.cdll.LoadLibrary(self.name_of_library)
                    self.potential_function = getattr(self.library, self.name_of_function)
                    os.chdir(self.work_dir)
                except AttributeError:
                    print("Didn't find the function in the provided shared library. Make sure the library path and potential function name are correct.")
                    sys.exit(1)
            else:
                print("We need either a library name or a module name. Please provide one.")
                sys.exit(1)

class TTM(Potential):
    def __init__(self, path_to_library: str, name_of_function="ttm_from_f2py", name_of_module="ttm", model=21):
        """Evaluates the energy and gradients of the TTM family of potentials.

        Args:
            model (int, optional): The TTM model which will be used. Options are 2, 21, and 3. Defaults to 21.
        """
        super().__init__(path_to_library=path_to_library, name_of_function=name_of_function, name_of_module=name_of_module)
        self.model = model
        self.initialize_potential()
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
        # Sadly, we need to re-order the geometry to TTM format which is all oxygens first.
        coords = self.ttm_ordering(coords)
        os.chdir(self.path_to_library)
        gradients, energy = self.potential_function(self.model, np.asarray(coords).T, int(len(coords) / 3))
        os.chdir(self.work_dir)
        return energy / 627.5, (-self.normal_water_ordering(gradients.T) / 627.5) / 1.88973
    
    def __call__(self, coords):
        return self.evaluate(coords)

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['potential_function']
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__dict__.update({"potential_function": None})
        self.initialize_potential()

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

class MBPol(Potential):
    def __init__(self, path_to_library: str, name_of_function="calcpotg_", name_of_library="./libmbpol.so"):
        super().__init__(path_to_library=path_to_library, name_of_function=name_of_function, name_of_library=name_of_library)

    def initialize_potential(self, num_waters):
        super().initialize_potential()
        if self.potential_function.argtypes is None:
            self.potential_function.restype = None
            self.potential_function.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
            ]
        self.num_waters = num_waters
        self.c_num_waters = ctypes.byref(ctypes.c_int32(self.num_waters))
    def evaluate(self, coords):
        if isinstance(coords, np.ndarray):
            self.initialize_potential(coords.shape[0] // 3)
        else:
            self.initialize_potential(len(coords) // 9) # N_w X 3 X 3
            coords = np.array(coords)
        coords=np.ascontiguousarray(coords, dtype=np.float64).flatten()
        grads = np.zeros_like(coords)
        potential_energy = np.zeros(1)
        self.potential_function(self.c_num_waters, potential_energy, coords, grads)
        return potential_energy[0] / 627.5, -np.reshape(grads, (3 * self.num_waters, 3)) / 627.5 / 1.88973
    
    def __call__(self, coords):
        return self.evaluate(coords)

class Protonated_Water(Potential):
    def __init__(self, num_waters: int, library_path: str, do_init=True):
        self.num_waters = num_waters
        self.library_path = library_path
        self.work_dir = os.getcwd()
        sys.path.insert(0, self.library_path)
        os.chdir(self.library_path)
        self.module = importlib.import_module("Protonated_Water")
        self.energy_function = getattr(self.module, "get_energy")
        self.energy_and_gradient_function = getattr(self.module, "get_energy_and_gradients")
        self.init_function = getattr(self.module, "initialize_potential")
        if do_init:
            self.init_function(self.num_waters)
        os.chdir(self.work_dir)

    def evaluate(self, coords, get_gradients=True):
        if get_gradients:
            return self.get_energy(coords)
        else:
            return self.get_energy_and_gradients(coords)

    def get_energy(self, coords):
        """
        Gets potential energy in hartree from coords.
        """
        os.chdir(self.library_path)
        energy = self.energy_function(coords.T)
        os.chdir(self.work_dir)
        return energy
    
    def get_energy_and_gradients(self, coords):
        """
        Gets potential energy in hartree and gradients in hartree per bohr from coords.
        """
        os.chdir(self.library_path)
        energy, gradients = self.energy_and_gradient_function(coords.T)
        os.chdir(self.work_dir)
        gradients = np.reshape(gradients, np.shape(coords))
        return energy, gradients