import xtb
from xtb import GFN0

import sys

import ase
from ase.io import read, write
from ase import Atoms
from ase.build import surface, add_adsorbate, molecule, add_vacuum
from ase.constraints import FixAtoms, FixedPlane
from ase.units import Hartree
from ase.constraints import ExpCellFilter
from ase.optimize.lbfgs import LBFGS
from ase.constraints import FixAtoms

inpName = sys.argv[1]
# read molecular structure data, here from a VASP geometry input
slab = read(inpName)

# create the calculator for GFN0-xTB under periodic boundary conditions
slab.calc = GFN0()

# initial single point calculation
e = slab.get_potential_energy()
print("Initial energy: eV, Eh", e, e/Hartree)

# setup optimization of cell parameters
#ecf = ExpCellFilter(mol)
#precon = Exp(A=3)
#relax = PreconFIRE(mol, precon=precon, trajectory='xtbopt.traj')

relax = LBFGS(slab,
maxstep=0.1)
#logfile='log',
#trajectory='ase.traj')
# do the optimization
relax.run(fmax=0.03)

# get the final single point energy
e = slab.get_potential_energy()
print("Final energy:   eV, Eh", e, e/Hartree)

# write final geometry to file
write('out-'+inpName, slab, format='vasp')
