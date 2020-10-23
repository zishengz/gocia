
import sys, os
import numpy as np
import ase.io as fio
from ase.db import connect
from gocia.utils.dbio import get_traj
import gocia.ensemble.comparator as comp
from gocia.ensemble import clusterIsomer

trajDBName = sys.argv[1]
geomCutoff = 0.1
enerCutoff = 0.1

print(' * Extracting data from database: %s'%trajDBName)
rawDB = connect(trajDBName)
traj = get_traj(rawDB.select())
ene = [img.info['eV'] for img in traj]
mag = [img.info['mag'] for img in traj]

# eneDiff = comp.compare_ene(ene, enerCutoff)
# geomSim = comp.compare_geom(traj, geomCutoff)
# bothPass = comp.bothSim(geomSim, eneDiff)

eneDiff = comp.compare_ene(ene, enerCutoff)
allEigDist = comp.compare_posEig(traj, geomCutoff)
bothPass = comp.bothSim(allEigDist, eneDiff)

# $ python test.py raw-N1.db
# # Timing: 3.050 s
# bothPass = comp.compare_dual(traj, geomCutoff, enerCutoff)
clusterIsomer(traj, bothPass, outName='sort-'+trajDBName.split('.')[0])

