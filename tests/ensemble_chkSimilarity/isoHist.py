import sys, os
import numpy as np
import ase.io as fio
from ase.db import connect
from gocia.utils.dbio import get_traj
from gocia.utils.visualize import histogram

trajDBName = sys.argv[1]
geomCutoff = 0.9
enerCutoff = 0.1

print(' * Extracting data from database: %s'%trajDBName)
rawDB = connect(trajDBName)
traj = get_traj(rawDB.select())
ene = np.array([img.info['eV'] for img in traj])

histogram(ene, traj[0].get_chemical_formula())