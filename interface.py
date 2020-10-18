#
# Created on Wed Oct 14 2020
#
# Copyright (c) 2020 Zisheng Zhang
# Alexandrova Lab, UCLA, United States

"""
This module defines the Surface object. 
"""

from gocia.data import elemSymbol, covalRadii
import ase.io as fio
import numpy as np
import json

class Interface:
    def __init__(self,
                 tags=None,
                 allAtoms=None,
                 subAtoms=None,
                 allPos=None,
                 constraint=None,
                 fixList=None,
                 bufferList=None,
                 cellParam=None,
                 pbcParam=None,
                 zLim=None,
                 adsList=None
                 ):

        interface = None

        if subAtoms is not None:
            if type(subAtoms) is str:
                self.subAtoms   = fio.read(subAtoms)
            else:
                self.subAtoms   =subAtoms
            self.constraints= self.subAtoms.constraints
            self.fixList    = self.constraints[0].get_indices()
            self.bufferList = [i for i in list(range(len(subAtoms)))\
                if i not in self.fixList]
            self.cellParam  = self.subAtoms.get_cell()
            self.pbcParam   = self.subAtoms.get_pbc()            
        
        if allAtoms is not None:
            if type(allAtoms) is str:
                self.allAtoms   = fio.read(allAtoms)
            else:
                self.allAtoms   = allAtoms
            self.allPos     = self.allAtoms.get_positions()

        if zLim is None:
            self.zLim       = [min(self.allPos[:,2][self.bufferList]) + 0.1,
                    max(self.allPos[:,2][self.bufferList]) + 2.5]

        if adsList is None:
            self.adsList    = np.array(
                [[self.allAtoms.get_atomic_numbers()[i], i]\
                for i in list(range(len(self.allAtoms)))\
                if i not in list(range(len(self.subAtoms)))]
            )

        if tags is not None:
            self.tags = tags
        else:
            self.tags = str(self.allAtoms)

        self.update()


    def __len__(self):
        return len(self.allAtoms)
    
    def copy(self):
        import copy
        myCopy = self.__class__(
            tags=self.tags,
            subAtoms=self.subAtoms,
            allAtoms=self.allAtoms
        )
        return myCopy

    def get_allAtoms(self):
        return self.allAtoms.copy()

    def set_allAtoms(self, newAllAtoms):
        newAllAtoms.wrap()
        self.allAtoms = newAllAtoms
        self.update()
    
    def wrap(self):
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.wrap()
        self.set_allAtoms(tmpAtoms)

    def get_cell(self):
        return self.cellParam.copy()

    def get_positions(self):
        return self.get_allAtoms().get_positions().copy()

    def set_positions(self, newPos):
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.set_positions(newPos)
        self.set_allAtoms(tmpAtoms)

    def get_chemical_symbols(self):
        return self.get_allAtoms().get_chemical_symbols().copy()

    def get_adsList(self):
        if len(self.get_allAtoms()) == len(self.subAtoms):
            return []
        else:
            return np.array(
                [[self.allAtoms.get_atomic_numbers()[i], i]\
                for i in list(range(len(self.allAtoms)))\
                if i not in list(range(len(self.subAtoms)))]
            )
    
    def get_adsIndexList(self):
        return list(self.get_adsList()[:,1])

    def get_topLayerList(self, depth = 1):
        allZ = self.allPos[:,2]
        return [a.index for a in self.allAtoms\
            if max(allZ) - allZ[a.index] <= depth]

    def get_bufferList(self):
        return self.bufferList.copy()

    def merge_adsorbate(self, adsAtoms):
        self.set_allAtoms(
            self.get_allAtoms().extend(adsAtoms)
        )

    def remove_adatom(self, adsIndexList=None):
        oldAdsList = self.get_adsIndexList()
        if adsIndexList is None:
            adsIndexList = list(range(len(oldAdsList)))
        if type(adsIndexList) is not list:
            adsIndexList = [adsIndexList]
        adsIndexList = [oldAdsList[i] for i in adsIndexList]
        tmpAtoms = self.get_allAtoms()
        print(len(tmpAtoms))        
        del tmpAtoms[adsIndexList]
        print(len(tmpAtoms))   
        self.set_allAtoms(tmpAtoms)

    def rattle(self, stdev = 0.1):
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.rattle(stdev = stdev, seed=np.random.seed())
        self.set_allAtoms(tmpAtoms)


    def preopt_lj(self, fileBaseName='tmp',\
        toler=0.2, stepsize=0.05, nsteps=200):
#        from ase.calculators.lj import LennardJones
        from gocia.calc.lj import LennardJones
        from ase.optimize.bfgs import BFGS
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.calc = LennardJones(tolerAngs=toler, tolerMult=toler)
        geomOpt = BFGS(
            tmpAtoms,
            maxstep=stepsize,
            trajectory=fileBaseName+'.traj',
            logfile=fileBaseName+'.log'
        )
        geomOpt.run(fmax = 0.01, steps = nsteps)
        geomOpt = fio.read(fileBaseName+'.traj', index=':')
        print('L-J pre-optimization converged in %i loops'%(len(geomOpt)))
        tmpAtoms = geomOpt[-1]
        self.set_allAtoms(tmpAtoms)

    def preopt_hooke(self, cutoff = 1.5,
        toler=0.2, stepsize=0.05, nsteps=200):
#        from ase.calculators.lj import LennardJones
        from gocia.calc.hooke import Hooke
        from ase.optimize.bfgs import BFGS
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.calc = Hooke(
            cutoff=cutoff,
            tolerAngs=toler,
            tolerMult=toler
            )
        geomOpt = BFGS(
            tmpAtoms,
            maxstep=stepsize,
            trajectory=None,
            logfile=None
        )
        geomOpt.run(fmax = 0.01, steps = nsteps)
#        print(tmpAtoms.get_potential_energy())
#        print('Hooke pre-optimization converged in %i loops'%(len(geomOpt)))
        self.set_allAtoms(tmpAtoms)




    def update(self):
        self.allPos  = self.allAtoms.get_positions()
        self.allAtoms.set_cell(self.cellParam)
        self.allAtoms.set_pbc(self.pbcParam)
        self.allAtoms.set_constraint(self.constraints)
        self.adsList = self.get_adsList()
        self.bufferList = [i for i in list(range(len(self.subAtoms)))\
                           if i not in self.fixList]
        

    def print_info(self):
        print(self.tags)
        print(self.allAtoms)
        print('Buffer region:\tZ = %.3f to %.3f'%\
            (self.zLim[0], self.zLim[1]))
        if self.get_adsList() == []:
            print('No adsorbates so far.')
        else:    
            print('Adsorbate:\t'+\
                str([elemSymbol[i] for i in self.adsList[:,0]]))

    def write(self, fileName):
        fio.write(fileName, self.get_allAtoms())






        


