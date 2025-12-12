import faulthandler
import signal

faulthandler.enable()
# Optional: Dump the traceback when you send a signal (e.g., SIGUSR1 on Unix or SIGINT for Ctrl+C)
def dump_traceback(signal_number, frame):
    faulthandler.dump_traceback()
signal.signal(signal.SIGUSR1, dump_traceback)

import ase, ase.io
from ase.db import connect
from ase.optimize import BFGS
from ase.neighborlist import NeighborList, natural_cutoffs
from gocia.ga.popGrandCanonPoly import PopulationGrandCanonicalPoly
from time import sleep
import os
from gocia.utils.ase import geomopt_iterate
import numpy as np
import input
from ase.calculators.lj import LennardJones

mycalc = LennardJones()

base = ase.io.read("all-POSCAR")
base.calc = mycalc

totConf = 1000
minConf = 5
subsPot = base.get_potential_energy()

pop = PopulationGrandCanonicalPoly(
    gadb='ini.db',
    substrate='sub-POSCAR',
    popSize = input.popSize,
    convergeCrit=input.popSize*10,
    subsPot = subsPot,
    chemPotDict = input.chemPotDict,
    zLim = input.zLim,
    )

try:
    list_file =os.listdir('.')
    list_job = [f for f in list_file if f[0] == 's' and '.' not in f]
    list_job.remove('sub-POSCAR')
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

while not pop.is_converged() or kidnum < minConf:
    if 'STOP' in os.listdir():
            print('STOP REQUESTED')
            exit()
    if kidnum > totConf:
            print('MAX # SAMPLE REACHED')
            exit()

    # GENERATE CHILD
    kidnum += 1
    kiddir = f's{str(kidnum).zfill(6)}'
    cwd = os.getcwd()
    try:
        os.mkdir(kiddir)
    except:
        pass
    os.chdir(kiddir)

    kid = None
    while kid is None:
        kid = pop.gen_offspring_box(
            mutRate=1,
            xyzLims=np.array([[0, 11.9902680382192521], [0, 13.1151602652946266], [11, 17]]),
            bondRejList = [['H', 'H'], ['O', 'O']],
            constrainTop=False,
            transVec=[[-1,1,-2,2,-4,4],[-1,1,-2,2,-4,4]],
            clusterAtoms=['Rh','Rh'],
            rattleOn=False,
            growOn=False,
            leachOn=False,
            moveOn=False,
            moveClusterOn=True,
            permuteOn=False,
            transOn=False
            )

    tmpatoms = kid.get_allAtoms()
    tmpatoms.calc = mycalc
    
    pop.add_aseResult(tmpatoms, workdir=kiddir, fn_frag='fragments')
    pop.natural_selection()

    os.chdir('..')

print(f'Done')
    
