from ase.db import connect
import ase.io as fio
import sys, os

dbName = sys.argv[1]
baseName = os.path.basename(dbName)
rootName = os.path.splitext(baseName)[0]
with connect(dbName) as myDB:
    for r in myDB.select():
        thisJobDir = '%s_%06d'%(rootName, r.id)
        os.makedirs('tmp/' + thisJobDir)
        fio.write('tmp/%s/%s.vasp'%(thisJobDir, thisJobDir), r.toatoms())
    

