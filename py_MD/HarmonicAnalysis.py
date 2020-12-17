#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
import subprocess as sub
import itertools as itt
import sys
import time
import uuid
from pathlib import Path
from .Fragments import Fragments
from .Potential import *
from .MBE_Potential import MBE_Potential
from .read_geometries import read_geoms

class Constants:
    """Helps me keep track of constants and conversions,
     originally written by Mark (b3m2a1)."""
    atomic_units = {
        "wavenumbers" : 4.55634e-6,
        "angstroms" : 1/0.529177,
        "amu" : 1.000000000000000000/6.02213670000e23/9.10938970000e-28   #1822.88839  g/mol -> a.u.
    }

    masses = {
        "H" : ( 1.00782503223, "amu"),
        "O" : (15.99491561957, "amu"),
        "D" : (2.0141017778,"amu"),
        "C" : (11.9999999958,"amu"),
        "N" : (14.003074,"amu")
    }
    @classmethod
    def convert(cls, val, unit, to_AU = True):
        vv = cls.atomic_units[unit]
        return (val * vv) if to_AU else (val / vv)

    @classmethod
    def mass(cls, atom, to_AU = True):
        m = cls.masses[atom]
        if to_AU:
            m = cls.convert(*m)
        return m

class HarmonicAnalysis:
    def __init__(self,eqGeom,atoms,potential,dx=1.0e-3,ofile=None):
        self.eqGeom = eqGeom
        self.atoms = atoms
        self.potential = potential
        self.dx = dx
        self.ofile=ofile
        self.nEls = 3 * len(self.atoms)

    def genStencil(self,dispTup,dim):
        cds = self.eqGeom
        dx = self.dx
        if dim == 1:
            """Generates 5-point 1D stencil for finite difference"""
            atm = dispTup[0] #atom of interest
            cd = dispTup[1] #x,y, or z
            stnShape = np.concatenate(([5], np.shape(cds)))
            stencilCds = np.zeros(stnShape)
            dxOrder = [-2*dx, -dx, 0, dx,2*dx]
            ct = 0
            for disp in dxOrder:
                dxAr = np.zeros(cds.shape)
                dxAr[atm, cd] += disp
                stencilCds[ct] = cds + dxAr
                ct+=1
            vs = self.potential(stencilCds) #Takes in an nxnAtomsx3 numpy array
            return vs
        elif dim  == 2:
            """Generates 9-point 2D stencil for finite difference"""
            atm1 = dispTup[0]
            cd1 = dispTup[1]
            atm2 = dispTup[2]
            cd2 = dispTup[3]
            stnShape = np.concatenate(([9], np.shape(cds)))
            stencilCds = np.zeros(stnShape)
            dxOrder = [-dx, 0, dx]
            ct = 0
            for disp in dxOrder:
                for disp2 in dxOrder:
                    dxAr = np.zeros(cds.shape)
                    dxAr[atm1, cd1] += disp
                    dxAr[atm2, cd2] += disp2
                    stencilCds[ct] = cds + dxAr
                    ct += 1
            vs = self.potential(stencilCds)
            return vs.reshape(3, 3)


    def finiteDiff(self,stencil,dim):
        dx = self.dx
        """Weights are well known, but taken from http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf"""
        if dim == 1: #second derivative
            #5 point stencil
            wts = np.array([
                [-1.0, 16.0, -30.0, 16.0, -1.0],
            ]) / (12.0*dx**2)
            der = stencil.dot(wts[0])
            return der
        elif dim == 2:
            # 9 point stencil
            wts = np.array([
                [1., 0., -1.],
                [1., 0., -1.]
            ]) / (dx * 2.)
            der = (stencil.dot(wts[0])).dot(wts[1])
            return der

    def genHess(self):
        nAtoms = len(self.atoms)
        hess = np.zeros((self.nEls,self.nEls))
        ###### Off Diagonals #############
        d=np.ravel_multi_index(np.array((np.triu_indices(nAtoms*3, 1))),[nAtoms*3,nAtoms*3]) #indexing nonsense
        atmInfo = np.array((np.unravel_index(d, (nAtoms, 3, nAtoms, 3)))).T #indexing nonsense (getting the upper indices of a 4D array for all my displacements)
        #atmInfo has information (cd1,atm1,cd2,atm2) which tells the stencil function which atom and dimension to displace
        offDiags = list(itt.combinations(np.arange(self.nEls),2)) #indices of upper triangle of hessian
        for k in range(len(offDiags)):
            stencil = self.genStencil(atmInfo[k],dim=2)
            hess[offDiags[k]] = self.finiteDiff(stencil,dim=2)
        hess = hess+hess.T #hessian is symmetric matrix

        ###### On Diagonals #############
        s = [np.arange(nAtoms),np.arange(3)]
        atmInfo = list(itt.product(*s))
        for onDiags in range(self.nEls):
            stencil1D = self.genStencil(atmInfo[onDiags],dim=1)
            hess[onDiags,onDiags]=self.finiteDiff(stencil1D,dim=1)
        return hess

    def diagonalize(self,hessian):
        masses = np.array([Constants.mass(a) for a in self.atoms])
        massesDup = np.repeat(masses, 3)
        hessianMW = hessian / np.sqrt(massesDup) / np.sqrt(massesDup[:, np.newaxis])
        freqs, normalModes = la.eigh(hessianMW)
        freqsCM = Constants.convert(np.sqrt(freqs), 'wavenumbers', to_AU=False)
        print(freqsCM)  # should be the same as freqsCM, should have 3 negative frequencies. and 3 others close to zero.
        xyz = ['_x','_y','_z']
        atmOut =[]
        for i in self.atoms:
            for j in xyz:
                atmOut.append(i+j)

        timestr = time.strftime("%Y%m%d-%H%M%S")
        if self.ofile is None:
            np.savetxt("frequencies" + timestr + ".txt", freqsCM)
        else:
            np.savetxt(str(self.ofile) + "_freq.txt", freqsCM)
        #np.savetxt("normalModes.txt",np.column_stack((atmOut,normalModes)),fmt='%s')

    def run(self):
        hessian = self.genHess()
        self.diagonalize(hessian)


def partridgePot(cds):
    """Calls executable calc_h2o_pot, which takes in the file hoh_coord.dat and saves the energies to hoh_pot.dat
    hoh_coord.dat has the form
    nGeoms
    hx hy hz
    hx hy hz
    ox oy oz
    h'x h'y h'z
    h'x h'y h'z
    o'x o'y o'z
    h''x h''y h''z
    h''x h''y h''z
    o''x o''y o''z
    ...
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Where nGeoms is an integer corresponding to the number of geometries you pass in
    """
    np.savetxt("PES0/hoh_coord.dat", cds.reshape(cds.shape[0] * cds.shape[1], cds.shape[2]), header=str(len(cds)),
               comments="")
    sub.run('./calc_h2o_pot', cwd='PES0')
    return np.loadtxt("PES0/hoh_pot.dat")-(  -1.9109019308531233E-006)

def TTM21FPot(cds):
    '''
    Calls executable TTM21F_energy.x, which takes in file hoh_coord.dat and saves the energies to hoh_pot.dat.
    The potential expects water molecules in O H H order.
    '''
    grid_energies = []
    #write hoh_coord.dat
    for iGeom, atoms__ in enumerate(cds):
        N=len(cds[iGeom])
        atoms__ = np.copy(cds[iGeom])
        atoms__ = [str(atoms__[x]*.529177) for x in range(len(atoms__))]
        for i in range(len(atoms__)):
            if i % 3 == 0:
                atoms__[i] = "O " + atoms__[i]
            else:
                atoms__[i] = "H " + atoms__[i]
        hoh_coord = str(uuid.uuid4())
        hoh_pot = str(uuid.uuid4())
        with open(hoh_coord, 'w') as f:
            f.write(str(N) + "\n")
            f.write("\n")
            for line in atoms__:
                printable = line.replace("[","").replace("]","")
                f.write(printable)
                f.write("\n")
        #evaluate potential and return energy
        with open(hoh_pot, 'w') as outfile:
            sub.call(["TTM21F_energy.x", hoh_coord], stdout=outfile)
            proc=sub.Popen(["grep","Energy", hoh_pot],stdout=sub.PIPE)
            energy = proc.stdout.read()
            grid_energies.append(float(energy.split()[3])*0.0015936011) #convert kcal/mol to a.u.
        sub.call(["rm", hoh_coord])
        sub.call(["rm", hoh_pot])
    return np.array(grid_energies)

class Potential_Wrapper:
    """
    Wraps an existing function that only calls a single geometry to returns many energies.
    """
    def __init__(self, pot_function, to_angstrom=True):
        self.pot_function = pot_function
        self.to_angstrom = to_angstrom
    
    def evaluate(self, cds):
        grid_energies = []
        for atoms in cds:
            grid_energies.append(self.pot_function(atoms / 1.88973)[0])
        return np.asarray(grid_energies)

if __name__ == '__main__':
    try:
        infile = sys.argv[1]
    except:
        print ("Expecting xyz formatted input file as command line argument. Exiting.")
        sys.exit(1)
    dxx = 1.e-2
    """Everything is in  Atomic Units going into generating the Hessian."""
    #read in all the geometries
    #header, labels, geoms = read_geoms(infile)
    fragments = Fragments(infile)
    geom = np.vstack(fragments.fragments)
    #ttm21f = TTM(["ttm*"], "ttm", "ttm_from_f2py", 21)
    mbpol = MBPol()
    mbe_ff = MBE_Potential(5, fragments, mbpol, return_extras=False)
    pot = Potential_Wrapper(mbe_ff.evaluate_on_geometry)
    # Perform harmonic analysis on geometry
    geom = Constants.convert(geom, \
    "angstroms",to_AU=True) #To Bohr from angstroms
    #atoms = labels[iGeom]
    atoms = ["O", "H", "H"] * 10
    HA_h2o = HarmonicAnalysis(eqGeom=geom,
                            atoms=atoms,
                            potential=pot.evaluate,
                            dx=dxx,
                            ofile=Path(infile)
                            )
    HarmonicAnalysis.run(HA_h2o)
