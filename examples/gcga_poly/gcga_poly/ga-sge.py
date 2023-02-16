
#from gocia.ga.popGrandCanon import PopulationGrandCanonical
from gocia.ga.popGrandCanonPoly import PopulationGrandCanonicalPoly
from gocia.hpc.sge import SGE
from time import sleep
import os
import datetime

import input

print(os.getpid())
nworker = 45
totConf = 2000
minConf = 1000

print('read in population')
pop = PopulationGrandCanonicalPoly(
    gadb='gcga.db',
    substrate='substrate.vasp',
    popSize = input.popSize,
    convergeCrit=input.popSize*10,
    subsPot = input.subsPot,
    chemPotDict = input.chemPotDict,
    zLim = input.zLim
    )

try:
    list_file =os.listdir('.')
    list_job = [f for f in list_file if f[0] == 's' and '.' not in f]
    list_jid = [int(j[1:]) for j in list_job]
    kidnum = max(list_jid)
    print(f'Restarting the search from kidnum = {kidnum}')
except:
    kidnum = 0
    print('New search! Starting from kidnum = 0')

if kidnum == 0:
    print('Initializing the population...')
    pop.initializeDB()
    pop.natural_selection()

print('Initializing the job queue...')
queue = SGE('zisheng')
while not pop.is_converged() or kidnum < minConf:
    while len(queue) < nworker:
        if 'STOP' in os.listdir():
            print('STOP REQUESTED')
            exit()
        if kidnum > totConf:
            print('MAX # SAMPLE REACHED')
            exit()
        kidnum += 1
        print('Job %i\t@%s'%\
            (kidnum, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        queue.submit('./sge-vasp.sh %06d'%kidnum)
        sleep(10*3)
    sleep(60*10)
    queue.update()
    queue.write()

if pop.is_converged():
    print('CONVERGED!')
else:
    print('TERMINATED!')
