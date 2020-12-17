from .Fragments import Fragments
from .Potential import *
from .MBE_Potential import MBE_Potential
from .Compute_Hessian import Hessian
import numpy as np
import sys, time

class Optimize:
    """Simple gradient descent implementation which works by taking a potential_function
    and following the gradient until certain convergence criteria are met.
    """
    def __init__(self, initial_geometry, 
                       potential_function, 
                       max_iterations=2000, 
                       max_force=10**-5, 
                       max_rms_force=10**-5, 
                       max_delta_energy=10**-9,
                       step_size=0.5):
        self.initial_geometry = initial_geometry
        self.potential_function = potential_function

        self.max_iterations = max_iterations

        # convergence conditions
        self.max_force = max_force
        self.max_rms_force = max_rms_force
        self.max_delta_energy = max_delta_energy

        # convergence parameters
        self.delta_energy = 10.0
        self.current_max_force = 10.0
        self.current_rms_force = 10.0

        self.step_size = step_size

    def hybrid_method(self):
        geometry = self.gradient_descent(stop_early=True)
        self.newtons_method(starting_geometry=geometry)

    def newtons_method(self, starting_geometry=None):
        if starting_geometry is None:
            geometry = np.copy(self.initial_geometry.flatten())
        else:
            geometry = starting_geometry.flatten()

        step_size = 1.0

        hessian_calculator = Hessian(self.potential_function)
        old_energy, gradients = self.potential_function(np.reshape(geometry, self.initial_geometry.shape))
        for iteration in range(self.max_iterations):
            hessian = hessian_calculator.evaluate(geometry)
            inverted_hessian = np.linalg.inv(hessian)
            geometry -= step_size * np.dot(inverted_hessian, gradients.flatten() / 1.88973)
            if abs(self.delta_energy) > self.max_delta_energy or self.current_max_force > self.max_force or self.current_rms_force > self.max_rms_force:
                energy, gradients = self.potential_function(np.reshape(geometry, self.initial_geometry.shape))

                self.update_convergence_parameters(energy, old_energy, gradients)
                old_energy = energy
                print(f"Iteration {iteration}: Energy: {energy*627.5:.6f}, ({self.delta_energy*627.5:.6f}); Max Force: {self.current_max_force:.6f}; RMS Force: {self.current_rms_force:.6f}")
            else:
                print("Converged Geometry:")
                geometry = np.reshape(geometry, self.initial_geometry.shape)
                print(geometry - np.mean(geometry, axis=0))
                return
    
    def gradient_descent(self, stop_early=False, stop_early_iteration=50):
        geometry = np.copy(self.initial_geometry)

        old_energy, f = self.potential_function(geometry)
        g = np.copy(f / 1.88973)
        h = np.copy(f / 1.88973)
        for iteration in range(self.max_iterations):
            old_geometry = np.copy(geometry)
            f_old = np.copy(f)
            f /= 1.88973
            geometry += self.step_size * f
            if abs(self.delta_energy) > self.max_delta_energy or self.current_max_force > self.max_force or self.current_rms_force > self.max_rms_force:
                energy, f_new = self.potential_function(geometry)
                if energy - old_energy > 0.0:
                    geometry = old_geometry
                    self.step_size *= 0.5
                    f = f_old
                    iteration -= 1
                    continue

                f_new /= 1.88973
                gamma = np.dot(f_new.flatten(), f_new.flatten()) / np.dot(g.flatten(), g.flatten())
                g = f_new
                h = g + gamma * h
                f = h

                self.update_convergence_parameters(energy, old_energy, f)

                old_energy = energy
                self.step_size *= 1.1
                print(f"Iteration {iteration}: Energy: {energy*627.5:.6f}, ({self.delta_energy*627.5:.6f}); Max Force: {self.current_max_force:.6f}; RMS Force: {self.current_rms_force:.6f}")
            else:
                print(f"Converged Geometry: Final Energy = {energy*627.5:.6f}")
                return geometry
            
            if stop_early is True and iteration >= stop_early_iteration:
                print(f"Stopped after {stop_early_iteration} iterations! Did not converge!")
                return geometry

        print(f"Failed to converge in {self.max_iterations} steps!")
    
    def update_convergence_parameters(self, current_energy, old_energy, gradients):
        self.delta_energy = current_energy - old_energy
        self.current_max_force = np.amax(np.abs(gradients))
        self.current_rms_force = np.sqrt(np.mean(np.einsum('ij,ij->i', gradients, gradients)))

if __name__ == '__main__':
    try:
        ifile = sys.argv[1]
    except:
        print("Didn't get an xyz file.")
        sys.exit(1)
    
    fragments = Fragments(ifile)
    #ttm21f = TTM(["ttm*"], "ttm", "ttm_from_f2py", 21)
    mbpol = MBPol()
    mbe_ff = MBE_Potential(5, fragments, mbpol, return_extras=False)

    optimizer = Optimize(np.vstack(fragments.fragments), mbe_ff.evaluate_on_geometry)

    start = time.time()
    print(optimizer.gradient_descent())
    #optimizer.newtons_method()
    #optimizer.hybrid_method()
    print(time.time() - start)