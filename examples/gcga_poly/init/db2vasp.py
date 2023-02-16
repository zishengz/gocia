
from ase.constraints import FixAtoms
from ase.io import read, write
import sys, os
import numpy as np

inp = sys.argv[1]

from ase.db import connect

with connect(inp) as db:
    for i in range(len(db)):
        curName = 's%s.vasp' % (str(i).zfill(6))
        curFrag = 's%s.fragments' % (str(i).zfill(6))
        struct = db.get(id=i+1).toatoms()
        write(curName, struct, vasp5=True)
        with open(curFrag, 'w') as f: 
            f.write(db.get(id=i+1).get('adsFrags'))
        print(' > %s written!' % curName)

