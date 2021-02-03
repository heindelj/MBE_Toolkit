from .Potential import Potential
from .MBE_Potential import MBE_Potential
from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms
from ase.units import Hartree, Bohr
import sys

class PotentialCalculator(Calculator):
    """
    Wraps the py_MD Potential object as an ASE calculator returing energy and forces
    """
    implemented_properties = ['forces', 'energy']
    nolabel = True

    def __init__(self, potential: Potential):
        super().__init__()
        self.pymd_potential = potential

    def calculate(self, atoms=None, properties=None, system_changes=all_changes,):
        if properties is None:
            properties = self.implemented_properties

        if type(atoms) is not Atoms:
            print("Can only evaluate energy of an ASE Atoms object.")
            sys.exit(1)
        # call the base class to initialize things
        Calculator.calculate(self, atoms, properties, system_changes)
        
        # Call pymd potential
        energy, gradients = self.pymd_potential.evaluate(atoms.get_positions())
        self.results['energy'] = energy * Hartree
        self.results['forces'] = -gradients * Hartree / Bohr
        return energy, -gradients

class MBEPotentialCalculator(Calculator):
    """
    Wraps the py_MD Potential object as an ASE calculator.
    """
    implemented_properties = ['forces', 'energy']
    nolabel = True

    def __init__(self, mbe_potential: MBE_Potential):
        super().__init__()
        self.mbe_potential = mbe_potential

    def calculate(self, atoms=None, properties=None, system_changes=all_changes,):
        if properties is None:
            properties = self.implemented_properties

        # call the base class
        Calculator.calculate(self, atoms, properties, system_changes)

        # Call pymd potential and convert to ase internal units
        output = self.mbe_potential.evaluate_on_geometry_parallel(atoms.get_positions())
        self.results['energy'] = output[0] * Hartree
        self.results['forces'] = output[1] * Hartree / Bohr
        return output[0], output[1]