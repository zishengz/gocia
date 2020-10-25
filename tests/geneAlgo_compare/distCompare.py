

import sys, os
import numpy as np
import ase.io as fio
from ase.db import connect
from gocia.utils.dbio import get_traj
import gocia.ensemble.comparator as comp
from gocia.ensemble import clusterIsomerAbs

trajDBName = sys.argv[1]
geomCutoff = 0.1
enerCutoff = 0.1

print(' * Extracting data from database: %s'%trajDBName)
rawDB = connect(trajDBName)
traj = get_traj(rawDB.select())
ene = [img.info['eV'] for img in traj]

# simMat = comp.compAll_srtDist_zz(traj)
# eneDiff = comp.compare_ene(ene, enerCutoff)
# bothPass = comp.bothSim(simMat, eneDiff)
# clusterIsomer(traj, bothPass, outName='sort-'+trajDBName.split('.')[0])

clusterIsomerAbs(traj)