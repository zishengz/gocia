import os
import numpy as np
from ase.db import connect
from ase.io import read, write

def get_keyVal(dbRows, keyVal):
    '''
    Returns 1-D list
    '''
    tmp = []
    for r in dbRows:
        tmp.append(r[keyVal])
    return tmp

def get_traj(dbRows):
    tmp = []
    for r in dbRows:
        a = r.toatoms()
        a.info = r.key_value_pairs
        tmp.append(a)
    return tmp

def get_projName():
    return os.getcwd().split('/')[-1]

def calypso2db():
    os.system('cd results;cak.py -n 9999 --vasp')
    eneData = [eval(r[-1]) for r in open('results/Analysis_Output.dat').readlines().split()]
    with connect('calypso-%s.db'%get_projName(), append=False) as myDb:
        for i in len(eneData):
            s = read('results/dir_origin/OCell_%i.vasp'%i+1)
            s.wrap()
            myDb.write(
                s,
                eV = eneData[i],
                done = 1,
            )
    print(' %i candidates writen to calypso-%s.db'%(len(eneData), get_projName()))
            
def vasp2db(nameKey=''):
    print(' --- Collecting vasp results ---')
    vdirs = [d.split('/')[0] for d in os.popen('ls *%s*/OSZICAR'%nameKey).readlines()]
    count_fin = 0
    with connect('vasp-%s.db'%get_projName(), append=False) as myDb:
        for d in vdirs:
            print('%s'%d, end='\t')
            info = open(d+'/OSZICAR', 'r').readlines()[-1].split()
            mag = eval(info[-1])
            ene_eV = eval(info[4])
            s = ai.read(d+'/CONTCAR')
            s.wrap()
			myDb.write(
				s,
				eV=ene_eV,
				mag=mag,
                done=1,
                )
            count_fin += 1
    print('\n %i candidates writen to vasp-%s.db'%(len(eneData), get_projName()))

def db2vasp(dbName):
    traj = read(dbName, index=':')
    for i in range(len(traj)):
        curName = 's%s.vasp'%(str(i).zfill(6))
        write(curName, traj[i], vasp5=True)
        print(' > %s written!'%curName)
