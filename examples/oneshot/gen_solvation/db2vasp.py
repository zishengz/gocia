
from ase.constraints import FixAtoms
from ase.io import read, write
import sys, os

inp = sys.argv[1]

from ase.db import connect

with connect(inp) as db:
    for i in range(len(db)):
        curName = 's%s.vasp' % (str(i).zfill(3))
        write(curName, db.get(id=i+1).toatoms(), vasp5=True)
        print(' > %s written!' % curName)

