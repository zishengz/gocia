from ase.db import connect
from gocia.utils.dbio import get_traj, get_projName
from gocia.ensemble import clusterIsomerAbs_GC
import matplotlib.pyplot as plt
import sys
from datetime import datetime

dbName = sys.argv[1]

with connect(dbName) as rawDB:
    traj = get_traj(rawDB.select())[:]
    grandPot = [img.info['grandPot'] for img in traj][:]

minima = []
avg = []
for i in range(len(grandPot)):
    minima.append(min(grandPot[:i+1]))
    if i < 25:
        nowPop = grandPot[:i+1]
    else:
        nowPop = grandPot[i-25:i]
    avg.append(sum(nowPop)/len(nowPop))

plt.figure(figsize=[9.6,4.8])
plt.scatter(list(range(len(grandPot))), grandPot, edgecolor='k', alpha=0.5)
plt.plot(minima, linewidth=1, alpha=0.9, label='current minimum')
plt.plot(avg, linestyle=':', label='current average')
plt.legend()
plt.grid(alpha=0.2)
plt.xlabel('# samples', fontsize='x-large')
plt.ylabel('GC Free Energy (eV)', fontsize='x-large')
plt.title(f'Proj:{get_projName()}, @{datetime.now()}')
plt.savefig('progress.png', dpi=500, bbox_inches = "tight", transparent=True)

clusterIsomerAbs_GC(
        traj,
        eneToler=0.01,
        geomToler1=5e-4,
        geomToler2=0.25,
        outName='sort-'+dbName.split('.')[0]
    )
