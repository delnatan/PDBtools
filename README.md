# PDBtools
Simple python tools to read/write PDB files

Example of use:

from PDBtools.PDB_ import PDB

# Load PDB file
lyz = PDB('6LYZ.pdb')

# Compute Scattering Intensity 
lyz.calc_profiles(qmax=0.4, dq=0.001)

# Interatomic distances, r and Probability Distribution, P(r) are then stored under
# lyz.r and lyz.Pr
# Scattering vector and intensities are under lyz.q and lyz.Iq
