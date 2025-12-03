
import random
import os, sys
from ase.io import read, write
from ase.db import connect
from ase.constraints import FixAtoms
from gocia.interface import Interface
from gocia.geom import build
from gocia.utils import conv_formula2list
import numpy as np

nSample = int(sys.argv[1])

surf = Interface(
    allAtoms = './Al2O3_Pt.vasp',
    subAtoms = './Al2O3.vasp',
    zLim=[7.5, 12.5]
)
surf.print()

gen_init = []
while len(gen_init) < nSample:
    print(f'\nGENERATING STRUCTURE {len(gen_init)+1}')
    newsurf = surf.copy()
    newsurf.rattle(0.2)

    newsurf = build.boxSample_adatom(
        newsurf,
        ['Pt']*9,
        xyzLims=np.array([[0.5, 6.5], [3.5, 9], [12, 18]]),
        toler_BLmin = -0.4,
        toler_CNmin = 1,
        toler_CNmax = 12,
        # bondRejList = [['N', 'N']],
        bondMustList= [['Pt', 'Pt']],
        doShuffle=True,
        rattle=False
    )

    if newsurf is not None:
        newsurf.preopt_hooke(
        cutoff = 1.0,
        toler = 0.1
        )
        s = newsurf.get_allAtoms()
#        del s.constraints
        gen_init.append(s)

print(f'WRITING STRUCTURES TO init.db')
write(f'init.db', gen_init)

