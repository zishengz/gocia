import sys

from ase.io import read, write
from ase.units import Hartree
from ase.optimize.lbfgs import LBFGS
from xtb.ase.calculator import XTB
from ase.db import connect

inpName = sys.argv[1]

slab = read(inpName)

slab.calc = XTB(method='GFN1-xTB', max_iterations=1000)
relax = LBFGS(
    slab,
    maxstep=0.1,
#    logfile=None,
#    trajectory='ase.traj'
    )
relax.run(fmax=0.05)

# get the final single point energy
e = slab.get_potential_energy()
print("Final energy:   eV, Eh", e, e/Hartree)

# write final geometry to file
write('out-'+inpName, slab)
