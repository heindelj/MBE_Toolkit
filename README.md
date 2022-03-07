**MBE_Toolkit**

The purpose of this package is to allow for the calculation of the many-body expansion (MBE) of energies and gradients for arbitrary potential surfaces and combinations of potential surfaces.

MBE_Toolkit allows you to calculate the MBE up to order N and even combine up to N potentials (the most possible) to calculate each of the N-body terms with a different potential. As an example, one may wish to calculate the 1- and 2-body terms at the CCSD(T) level, but only use HF for the 3-body term. Then, one can use a classical potential to get the 4- through N-body terms. When calculating these 4- to N-body terms, MBE_Toolkit calculates this by taking the difference of the MBE up to 3-body and the full calculation on the composite system.
