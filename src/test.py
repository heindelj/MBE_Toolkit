from Fragments import Fragments
from read_geometries import read_geoms
from Potential import get_ASE_NWChem_Potential
from MBE_Potential import ASE_MBE_Potential

nwchem = get_ASE_NWChem_Potential('scf', 'sto-3g')
w20_frags = Fragments("data/W20_global_minimum_ttm21f_fragmented.xyz", nwchem)
mbe_nwchem = ASE_MBE_Potential(3, w20_frags)
print(mbe_nwchem.evaluate_on_fragments())
