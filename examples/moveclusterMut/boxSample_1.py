
import random
import os, sys
from ase.io import read, write
from ase.db import connect
from ase.constraints import FixAtoms
from gocia_new.interface import Interface
from gocia_new.geom import build
from gocia_new.utils import conv_formula2list
import numpy as np
from ase.calculators.lj import LennardJones

# Usage:
# $ python boxSample.py xxx.xyz 50 
# If too slow:
# $ nohup python -u boxSample.py xxx.xyz 50 &

mycalc = LennardJones()

nSample = int(sys.argv[1])

surf = Interface(
    tags = 'surf',
    allAtoms = 'all-POSCAR',
    subAtoms = 'sub-POSCAR',
    zLim=[10, 18]
)
surf.print()

adslist=['H','O','HO','CO']
atomlist=['H','O','O','C']

with connect('ini.db') as db:
    while len(db) < nSample:
        mutlist = []
        bondMustList = []
        print(f'\nGENERATING STRUCTURE {len(db)+1}')
        tmpsurf = surf.copy()
        nMut = random.randint(1, 2)
        for i in range(nMut):
            j = random.randint(0,3)
            mutlist.append(adslist[j])
        temp = ['Rh']
        temp.append(atomlist[j])
        print(mutlist)
        print(bondMustList)
        tmpsurf = build.boxSample_frag(
            tmpsurf,
            mutlist,
            xyzLims=surf.get_sampling_box(),
            toler_BLmax=1.0,
            toler_BLmin=0.4,
            toler_CNmax=3,
            toler_CNmin=0,
            doShuffle=True,
            constrainTop=False,
            rattle=False, rattleStdev=0.05, rattleZEnhance=False,
            ljopt=False, ljstepsize=0.01, ljnsteps=400,
            bondCutoff = 0.80
            )
        
        tmpatoms = tmpsurf.get_allAtoms()
        tmpatoms.calc = mycalc

        db.write(
                tmpsurf.get_allAtoms(), 
                done=1, 
                eV = tmpatoms.get_potential_energy(),
                adsFrags="{}".format(tmpsurf.get_fragList()))

print(f'WRITING STRUCTURES TO test.db')

