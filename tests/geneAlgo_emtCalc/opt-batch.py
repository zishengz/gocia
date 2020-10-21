import sys
from ase.db import connect
from gocia.utils.dbio import get_traj

from ase.calculators.emt import EMT
from ase.optimize.bfgs import BFGS

def opt_emt(atoms):
    tmpAtoms = atoms.copy()
    tmpAtoms.calc = EMT()
    geomOpt = BFGS(tmpAtoms, maxstep=0.05,trajectory=None,logfile=None)
    geomOpt.run(fmax=0.001, steps=200)
    return tmpAtoms

inpName = sys.argv[1]
oldDB = connect(inpName)
traj = get_traj(oldDB.select())

newDB = connect('out-'+inpName, append=False)
for i in range(len(traj)):
    print('Optimizing\t%i/%i'%(i, len(traj)))
    newDB.write(opt_emt(traj[i]))

