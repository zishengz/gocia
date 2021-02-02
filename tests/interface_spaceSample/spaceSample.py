from gocia.geom import rand_point_box
from gocia.geom.build import boxSample_adatom
from gocia.interface import Interface
from ase.io import read, write
from gocia.utils import conv_formula2list
import sys
import numpy as np
from ase.db import connect

fn_substrate = sys.argv[1]
adatoms = sys.argv[2]
nsamples = int(sys.argv[3])

s = read('./armchair.vasp')
mycell = np.array(s.get_cell(complete=True))
myLim = np.array([
    [0, mycell[0][0]],
    [0, mycell[0][0]],
    [7.0, 11.0]
    ])

# 1 M loops: ~4.2 s
# for i in range(1000000):
#     myPnt = rand_point_box(myLim)

surf = Interface(s, s, zLim = myLim[2])
surf.print()

with connect('%s-%s.db'%(fn_substrate.split('.')[0], adatoms), append=False) as db:
    for i in range(nsamples):
        print(i,end='\t')
        genAds = boxSample_adatom(
            surf, conv_formula2list(adatoms), myLim,
            rattle=True, rattleZEnhance=True)
        db.write(genAds.get_allAtoms())


