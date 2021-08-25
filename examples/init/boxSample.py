
import random
import os, sys
import ase.io as ai
from ase.db import connect
from gocia.interface import Interface
from gocia.geom import build
from gocia.utils import conv_formula2list
import numpy as np

# Usage:
# $ python boxSample.py substrate.vasp 50 
# If too slow:
# $ nohup python -u boxSample.py substrate.vasp 50 &
surfName = sys.argv[1]
nSample = int(sys.argv[-1])

surf = Interface(
    tags = surfName.split('.')[0],
    allAtoms = ai.read(surfName),
    subAtoms = ai.read(surfName),
    zLim=[14, 21]
)
surf.print()

myDB = connect('tmp.db', append=False)
for i in range(nSample):
    # Set here the range of coverage
    nB = random.randint(16,20)
    nO = random.randint(int(nB*0.5),int(nB*1.5))
    newsurf = build.boxSample_adatom(
        surf,
        ['B']*nB+['O']*nO,
        xyzLims=np.array([[0,9.36],[0,9.36],[14,21]]),
        toler_BLmin = -0.2,
        toler_CNmax = 3,
        bondRejList = [['O', 'O']],
        doShuffle=True,
        rattle=True, rattleStdev=0.05,
    )
    newsurf.preopt_hooke(
       cutoff = 1.2,
       toler = 0.1
        )
    myDB.write(newsurf.get_allAtoms(), done=0)

os.system('mv %s %s'%('tmp.db', str(newsurf.get_allAtoms().symbols)+'.db'))

