
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

# # get the final single point energy
#e = slab.get_potential_energy()
#print("Final energy:   eV, Eh", e, e/Hartree)

# write final geometry to file
with connect('../done-%s.db'%inpName.split('_')[0]) as doneDB:
    doneDB.write(slab)

# Update into database
with connect('../../%s.db'%inpName.split('_')[0]) as myDB:
    myID = eval(inpName.split('.')[0].split('_')[1].lstrip('0'))
    myDB.write(slab, id=myID, eV=slab.get_potential_energy(), done=1)
