from ase.db import connect
from gocia.utils.dbio import get_traj, vasp2db, get_projName
from gocia.ensemble import clusterIsomerAbs

outName = 'vasp-%s.db'%get_projName()
vasp2db()

with connect(outName) as rawDB:
    traj = get_traj(rawDB.select())
    clusterIsomerAbs(
        traj,
        eneToler=0.05,
        geomToler1=1e-3,
        geomToler2=0.5,
        outName='sort-'+outName.split('.')[0]
    )
