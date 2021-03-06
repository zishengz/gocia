
import os, sys
import ase.io as ai
from ase.db import connect
from gocia.interface import Interface
from gocia.geom import build

# $ python directSample.py GaN_4-l_vac-1.vasp 'N N N' 10
# # Timing: 0.7767 s/structure
surfName = sys.argv[1]
adsInp = sys.argv[2].split()
nSample = int(sys.argv[3])

surf = Interface(
    tags = surfName.split('.')[0]+' + '+sys.argv[2],
    allAtoms = ai.read(surfName),
    subAtoms = ai.read(surfName),
    zLim=[11, 14]
)
surf.print()

myDB = connect('tmp.db', append=False)
for i in range(nSample):
    print('Structure %i'%i,end='\t')
    newsurf = build.grow_adatom(
        surf,
        adsInp,
        toler=0.25,
        sampZEnhance=4.0,
        doShuffle=True,
        cnCount=True,
        rattle=True, rattleStdev=0.1,
		rattleZEnhance=True,
        zLim=surf.zLim,
    )
    # newsurf.preopt_hooke(
    #    cutoff = 1.25,
    #    toler = 0.25
    #     )
    myDB.write(newsurf.get_allAtoms())
os.system('mv %s %s'%('tmp.db', str(newsurf.get_allAtoms().symbols)+'.db'))

