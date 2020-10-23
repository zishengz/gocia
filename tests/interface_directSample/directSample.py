
import os, sys
import ase.io as ai
from ase.db import connect
from gocia.interface import Interface
from gocia.geom import build
from gocia.utils import conv_formula2list

# $ python directSample.py GaN_4-l_vac-1.vasp 'N N N' 10
# # Timing: 0.7767 s/structure
surfName = sys.argv[1]
adsInp = conv_formula2list(sys.argv[2])
nSample = int(sys.argv[3])

surf = Interface(
    tags = surfName.split('.')[0]+' + '+sys.argv[2],
    allAtoms = ai.read(surfName),
    subAtoms = ai.read(surfName),
    zLim=None
)
surf.print()

myDB = connect('tmp.db', append=False)
for i in range(nSample):
    print('Structure %i'%i,end='\t')
    newsurf = build.grow_adatom(
        surf,
        adsInp,
        toler=0.25,
        doShuffle=True,
        cnCount=True,
        sameElemPenalty=0,
        rattle=True, rattleStdev=0.1,
		rattleZEnhance=True,
        zLim=surf.zLim
    )
    newsurf.preopt_hooke(
       cutoff = 1.2,
       toler = 0.1
        )
    myDB.write(newsurf.get_allAtoms(), done=0)
os.system('mv %s %s'%('tmp.db', str(newsurf.get_allAtoms().symbols)+'.db'))

