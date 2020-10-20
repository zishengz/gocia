
import sys, os
import numpy as np
import ase.io as fio
from ase.db import connect
from gocia.utils.dbio import get_traj
import gocia.ensemble.comparator as comp
from gocia.ensemble import clusterIsomer

trajDBName = sys.argv[1]
geomCutoff = 0.9
enerCutoff = 0.1

print(' * Extracting data from database: %s'%trajDBName)
rawDB = connect(trajDBName)
traj = get_traj(rawDB.select())
ene = [img.info['eV'] for img in traj]
mag = [img.info['mag'] for img in traj]

# # stepwise, with plot
# geomSim = detect.compare_geom(traj, geomCutoff)
# eneDiff = detect.compare_ene(ene, enerCutoff)
# bothPass = detect.bothSim(geomSim, eneDiff)
# from gocia.utils.plotter import heatmap
# heatmap(geomSim, 'geom'+str(geomCutoff))
# heatmap(eneDiff, 'ener'+str(enerCutoff))
# heatmap(bothPass, 'both_G%sE%s'%\
#     (str(geomCutoff), str(enerCutoff)))

# $ python test.py raw-N1.db
# # Timing: 3.050 s
bothPass = comp.compare_dual(traj, geomCutoff, enerCutoff)
clusterIsomer(traj, bothPass, outName='sort-'+trajDBName.split('.')[0])

