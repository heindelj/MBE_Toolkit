import sys, os, subprocess, importlib
import numpy as np
import ctypes

__all__ = ['Potential', 'TTM', 'MBPol']

class Potential:
    """Abstract base class for potential energy surfaces. A single function, evaluate(),
    must be implemented which returns the energy and gradients.

    Each child class of this needs to implement any methods to:
    1) returns the gradients with the same ordering as it received them
    2) help with evaluation of the potential
    """
    def __init__(self, dll_filenames=None, name_of_module=None, name_of_function=None):
        self.dll_filenames = dll_filenames
        self.name_of_module = name_of_module
        self.name_of_function = name_of_function
        self.work_dir = os.getcwd()
        self.pot_dir = os.path.dirname(__file__)
    
    def evaluate(self, coords):
        raise NotImplementedError

    def initialize_potential(self):
        """
            Imports a potential which is accessed via one or multiple dlls. The dlls should be stored in pyMD/bin.
            dll_filenames can be something like "ttm*" where "*" is a wildcard for each dll file.
            Also takes the name of python module, as a string, which calls the potential.

            This function is crucial if you want to use a multiprocessing pool, as it needs to be
            used as the initializer for the pool.

            dll_filenames: a list of filenames to all needed dlls or a keyword to get multiple of them.
            name_of_module: a string which is the name of a module that wraps the potential to be called.
            name_of_fucntion: a string which is the name of the function to be called from the imported module.
        """
        lib_path = os.path.join(self.pot_dir, "..", "bin")
        if type(self.dll_filenames) is not list:
            print("dll_filenames should be a list of files to access.")
            sys.exit(1)
        for file_name in self.dll_filenames:
            copy_command = "cp " + lib_path + os.path.sep + file_name + " " + self.pot_dir
            subprocess.call(copy_command, shell=True)
        os.chdir(self.pot_dir)
        os.environ['LD_LIBRARY_PATH'] = self.pot_dir
        try:
            module = importlib.import_module(self.name_of_module)
            pot_function = getattr(module, self.name_of_function)
            os.chdir(self.work_dir)
            return pot_function
            #for file_name in dll_filenames:
            #    remove_command = "rm " + os.path.join(self.pot_dir, file_name)
            #    subprocess.call(remove_command, shell=True)
        except ImportError:
            print("Did not find potential module. Make sure you have compiled it and the library can be linked against, including things like libgfortran and libgcc.")
            sys.exit(1)

class TTM(Potential):
    def __init__(self, dll_filenames: list, name_of_module: str, name_of_function: str, model=21):
        """Evaluates the energy and gradients of the TTM family of potentials.

        Args:
            model (int, optional): The TTM model which will be used. Options are 2, 21, and 3. Defaults to 21.
        """
        super().__init__(dll_filenames=dll_filenames, name_of_module=name_of_module, name_of_function=name_of_function)
        self.model = model
        self.pot_function = self.initialize_potential()
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
        os.chdir(self.pot_dir)
        gradients, energy = self.pot_function(self.model, np.asarray(coords).T, int(len(coords) / 3))
        os.chdir(self.work_dir)
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

class MBPol(Potential):
    lib_loc = os.path.join(os.path.dirname(__file__), "..","bin")
    def __init__(self):
        self._nw = None
        self._cnw = None
        self._v = np.zeros(1)
        self._lib = None
        self._fun = None
    def initialize(self, nw):
        if self._lib is None:
            cur_dir = os.getcwd()
            try:
                os.chdir(self.lib_loc)
                self._lib = ctypes.cdll.LoadLibrary("./libmbpol.so")
                self._fun = self._lib.calcpotg_
            finally:
                os.chdir(cur_dir)
            self._fun.restype = None
            self._fun.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
            ]
        self._nw = nw
        self._cnw = ctypes.byref(ctypes.c_int32(self._nw))
    def evaluate(self, coords):
        if isinstance(coords, np.ndarray):
            self.initialize(coords.shape[0] // 3)
        else:
            nw = len(coords) // 9 # N_w X 3 X 3
            self.initialize(nw)
            coords = np.array(coords)
        coords=np.ascontiguousarray(coords, dtype=np.float64).flatten()
        grads = np.zeros_like(coords)
        self._fun(self._cnw, self._v, coords, grads)
        return self._v[0] / 627.5, -np.reshape(grads, (3 * self._nw, 3)) / 627.5 / 1.88973
    def __call__(self, coords):
        return self.evaluate(coords)