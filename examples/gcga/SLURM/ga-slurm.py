
from gocia.ga.popGrandCanon import PopulationGrandCanonical
from gocia.hpc.slurm import SLURM
from time import sleep
import os
import datetime

import input

print(os.getpid())
nworker = 25
totConf = 2000
minConf = 1000

print('read in population')
pop = PopulationGrandCanonical(
    gadb='gcga.db',
    substrate='substrate.vasp',
    popSize = input.popSize,
    convergeCrit=input.popSize*10,
    subsPot = input.subsPot,
    chemPotDict = input.chemPotDict,
    zLim = input.zLim
    )

print('initializing queue')
queue = SLURM('zisheng')

# Comment if restart
pop.initializeDB()
pop.natural_selection()

kidnum = 0
while kidnum < totConf and\
    not pop.is_converged2()\
    and 'STOP' not in os.listdir()\
    or kidnum < minConf:
    while len(queue) < nworker:
        if 'STOP' in os.listdir(): exit()
        kidnum += 1
        print('Job %i\t@%s'%\
            (kidnum, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        queue.submit('./sge-vasp.sh %06d'%kidnum)
        sleep(10)
    sleep(300)
    queue.update()
    queue.write()

if pop.is_converged():
    print('CONVERGED!')
else:
    print('TERMINATED!')
