
from gocia.ga.popGrandCanon import PopulationGrandCanonical
from gocia.hpc.sge import SGE
from time import sleep
import os
import datetime

import input

print(os.getpid())
nworker = 30
totConf = 3000
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
queue = SGE('zisheng')

pop.initializeDB()
pop.natural_selection()

kidnum = 0
while kidnum < totConf and not pop.is_converged() and 'STOP' not in os.listdir()\
    or kidnum < minConf:
    while len(queue) < nworker:
        kidnum += 1
        print('Job %i\t@%s'%\
            (kidnum, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        queue.submit('./sge-vasp.sh %06d'%kidnum)
        sleep(10)
    sleep(60)
    queue.update()
    queue.write()

if pop.is_converged():
    print('CONVERGED!')
else:
    print('TERMINATED!')
