import sys
from ase.db import connect
from gocia.utils.dbio import get_traj

# fmax = 0.01, maxstep = 0.05  ==> 37.064 s
# fmax = 0.001, maxstep = 0.05 ==> 48.045 s
# fmax = 0.01, maxstep = 0.1   ==> 32.981 s
def opt_emt_bfgs(atoms):
    from ase.calculators.emt import EMT
    from ase.optimize.bfgs import BFGS
    tmpAtoms = atoms.copy()
    tmpAtoms.calc = EMT()
    geomOpt = BFGS(tmpAtoms, maxstep=0.1,trajectory=None,logfile=None)
    geomOpt.run(fmax=0.01, steps=1000)
    return tmpAtoms

# fmax = 0.01, maxstep = 0.1   ==> 28.852 s
def opt_emt_lbfgs(atoms):
    from ase.calculators.emt import EMT
    from ase.optimize.lbfgs import LBFGS
    tmpAtoms = atoms.copy()
    tmpAtoms.calc = EMT()
    geomOpt = LBFGS(tmpAtoms, maxstep=0.1,trajectory=None,logfile=None)
    geomOpt.run(fmax=0.01, steps=1000)
    return tmpAtoms


inpName = sys.argv[1]
oldDB = connect(inpName)
traj = get_traj(oldDB.select())

newDB = connect('out-'+inpName, append=False)
for i in range(len(traj)):
    print('Optimizing\t%i/%i'%(i+1, len(traj)))
    opt = opt_emt_lbfgs(traj[i])
    newDB.write(opt, eV=opt.get_potential_energy(), mag=0)

