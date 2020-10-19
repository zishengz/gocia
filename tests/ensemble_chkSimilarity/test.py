
import sys, os
import numpy as np
import ase.io as fio
from ase.db import connect
import gocia.geom.fingerprint as fgp
import gocia.utils.linalg as la
from gocia.utils.visualize import heatmap
from gocia.utils.dbio import get_keyVal, get_fgp
from gocia.utils.plotter import heatmap
from gocia.ensemble import clusterIsomer

# $ python test.py raw-N1.db
# # Timing: 6.042 s
trajDBName = sys.argv[1]
geomCutoff = 0.9
enerCutoff = 0.1

print(' * Extracting data from database: %s'%trajDBName)
rawDB = connect(trajDBName)
fgpArr = get_fgp(rawDB, fgp.coordinationFGPv1)
ene = np.array(get_keyVal(rawDB, 'eV'))
mag = np.array(get_keyVal(rawDB, 'mag'))
traj = fio.read(trajDBName, index=':')

print(' * Analyzing geometry similarity ...')
geomSim = la.cosSimMatrix(fgpArr)
geomSim[geomSim <  geomCutoff] = 0
geomSim[geomSim >= geomCutoff] = 1.0
heatmap(geomSim, 'geom'+str(geomCutoff))

print(' * Analyzing   energy similarity ...')
eneDiff = np.abs(la.diffMatrix(ene))
eneDiff[eneDiff <= enerCutoff] = -1.0
eneDiff[eneDiff >  enerCutoff] = 0
eneDiff[eneDiff == -1.0] = 1
heatmap(eneDiff, 'ener'+str(enerCutoff))

print(' * Analyzing  overall similarity ...')
bothPass = la.normalize_mat(eneDiff + geomSim)
heatmap(bothPass, 'both_G%sE%s'%\
    (str(geomCutoff), str(enerCutoff)))

print(' * Detecting unique isomers ...')
clusterIsomer(
    traj,
    bothPass,
    ene,
    mag,
)


