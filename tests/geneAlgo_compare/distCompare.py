

import sys
from ase.db import connect
from gocia.utils.dbio import get_traj
from gocia.ensemble import clusterIsomerAbs

trajDBName = sys.argv[1]

print(' * Extracting data from database: %s'%trajDBName)
with connect(trajDBName) as rawDB:
    traj = get_traj(rawDB.select())

clusterIsomerAbs(
    traj,
    eneToler=0.05,
    geomToler1=5e-4,
    geomToler2=0.5,
    outName='sort'
    )