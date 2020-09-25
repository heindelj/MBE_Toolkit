import numpy as np

class Hessian:
    """Computes the Hessian numerically if given a geometry and potential function.
    The potential function should return the energy only.
    """
    def __init__(self, potential_function):
        self.potential_function = potential_function

    def potential_wrapper(self, geometry, shape):
        return self.potential_function(geometry.reshape(shape))

    def evaluate(self, geometry, epsilon=1e-6):
        geometry = geometry.flatten()
        N = geometry.size
        hessian = np.zeros((N,N))
        df_0 = self.potential_wrapper(geometry, (-1, 3))[1].flatten()
        for i in range(N):
            xx0 = geometry[i]
            geometry[i] = xx0 + epsilon
            df_1 = self.potential_wrapper(geometry, (-1, 3))[1].flatten()
            hessian[i,:] = (df_1 - df_0)/epsilon
            geometry[i] = xx0
        return hessian

if __name__ == '__main__':
    from Fragments import Fragments
    from Potential import *
    from MBE_Potential import MBE_Potential
    import sys

    try:
        ifile = sys.argv[1]
    except:
        print("Didn't get an xyz file.")
        sys.exit(1)
    
    fragments = Fragments(ifile)
    ttm21f = TTM(21)
    mbe_ff = MBE_Potential(6, fragments, ttm21f)

    geometry = np.vstack(fragments.fragments)
    hessian_calculator = Hessian(mbe_ff.evaluate_on_geometry)
