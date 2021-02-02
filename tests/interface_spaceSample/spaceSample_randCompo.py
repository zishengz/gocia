from gocia.geom.build import boxSample_adatom
from gocia.interface import Interface
from ase.io import read, write
from gocia.utils import conv_formula2list
import sys
import numpy as np
from ase.db import connect
import random

fn_substrate = sys.argv[1]
nsamples = int(sys.argv[2])

s = read(fn_substrate)
mycell = np.array(s.get_cell(complete=True))
myLim = np.array([
    [0, mycell[0][0]],
    [0, mycell[0][0]],
    [7.0, 11.0]
    ])

surf = Interface(s, s, zLim = myLim[2])
surf.print()

with connect('%s-rand.db'%(fn_substrate.split('.')[0]), append=False) as db:
    for i in range(nsamples):
        nB = random.randint(2,8)
        nO = random.randint(int(nB/2),nB)
        nH = random.randint(1,int(nO/2)+1)
        adatoms = ['B']*nB + ['O']*nO + ['H']*nH
        print('%i\t%s'%(i, ''.join(adatoms)),end='\t')
        genAds = boxSample_adatom(
            surf, adatoms, myLim,
            toler_BLmin = -0.2,
            toler_BLmax = 0.2,
            toler_CNmin = 1,
            toler_CNmax = 3,
            rattle=True, rattleZEnhance=True)
        genAds.preopt_hooke(
        cutoff = 1.2,
        toler = 0.1
        )
        db.write(genAds.get_allAtoms(), done=0)


