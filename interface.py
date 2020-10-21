#
# Created on Wed Oct 14 2020
#
# Copyright (c) 2020 Zisheng Zhang
# Alexandrova Lab, UCLA, United States

"""
This module defines the Surface object. 
"""

from gocia.data import elemSymbol, covalRadii
from gocia import geom
import ase.io as fio
from ase.atoms import Atoms
from ase.constraints import FixAtoms
from ase.build.tools import sort
import numpy as np
import json

class Interface:
    def __init__(self,
                 allAtoms=None,
                 subAtoms=None,
                 fixList=None,
                 bufList=None,
                 adsList=None,
                 cellParam=None,
                 pbcParam=None,
                 zLim=None,
                 tags=None,
                 info=None
                 ):

        interface = None

        if subAtoms is not None:
            if type(subAtoms) is str:
                self.subAtoms = fio.read(subAtoms)
            else:
                self.subAtoms =subAtoms

        if allAtoms is not None:
            if type(allAtoms) is str:
                self.allAtoms   = fio.read(allAtoms)
            else:
                self.allAtoms   = allAtoms

        self.fixList    = self.subAtoms.constraints[0].get_indices()
        self.bufList = [i for i in list(range(len(subAtoms)))\
            if i not in self.fixList]
        self.adsList    = [i for i in list(range(len(self.allAtoms)))\
            if i not in list(range(len(self.subAtoms)))]

        self.cellParam  = self.subAtoms.get_cell()
        self.pbcParam   = self.subAtoms.get_pbc()            
        
        if zLim is None:
            allPos = self.allAtoms.get_positions()
            self.zLim = [min(allPos[:,2][self.bufList]) + 0.1,
                    max(allPos[:,2][self.bufList]) + 2.5]
        else:
            self.zLim = zLim

        if tags is not None:
            self.tags = tags
        else:
            self.tags = self.get_formula()

        if info is not None:
            self.info = info
        else:
            self.info={}

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

    def update(self):
        '''
        Call after updating the allAtoms!
        '''
        # Get buffer indices by removing fixed atoms from substrate
        self.bufList = [i for i in list(range(len(self.subAtoms)))\
                        if i not in self.fixList]
        # Get adsorbate indices by removing subAtoms from the allAtoms
        self.adsList = []
        if len(self.get_allAtoms()) != len(self.subAtoms):
            self.adsList = [i for i in list(range(len(self.allAtoms)))\
                            if i not in list(range(len(self.subAtoms)))]
        self.fixList = self.subAtoms.constraints[0].get_indices()
        self.allAtoms.set_constraint(FixAtoms(self.fixList))
        self.allAtoms.set_cell(self.get_cell())
        self.allAtoms.set_pbc(self.get_pbc())

    def sort(self):
        ads = self.get_adsAtoms().copy()
        ads = sort(ads)
        self.set_adsAtoms(ads)

    def print(self):
        print('#TAG: ', self.tags)
        print(' |-Fixed atoms:    ', self.get_fixAtoms().symbols)
        print(' |-Buffer atoms:   ', self.get_bufAtoms().symbols)
        if self.get_adsList() == []:
            print(' |-Adsorbates:     None')
        else:    
            print(' |-Adsorbates:     ',self.get_adsAtoms().symbols)
        print(' |-Buffer region:   Z = %.3f to %.3f'%\
            (self.zLim[0], self.zLim[1]))
        print(' |-Info:           ', self.info, '\n')
        
    def get_formula(self):
        return self.allAtoms.get_chemical_formula()

    def get_cell(self):
        return self.cellParam.copy()

    def get_pbc(self):
        return self.pbcParam.copy()

    def get_pos(self):
        return self.get_allAtoms().get_positions().copy()

    def get_subPos(self):
        return self.get_subAtoms().get_positions().copy()

    def get_fixList(self):
        return self.fixList.copy()

    def get_bufList(self):
        return self.bufList.copy()

    def get_adsList(self):
        return self.adsList.copy()

    def get_optList(self):
        return self.bufList.copy() + self.adsList.copy()

    def get_allAtoms(self):
        return self.allAtoms.copy()

    def get_subAtoms(self):
        return self.subAtoms.copy()

    def get_fixAtoms(self):
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index not in self.get_fixList()]]
        return tmpAtoms

    def get_bufAtoms(self):
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index not in self.get_bufList()]]
        return tmpAtoms
    
    def get_adsAtoms(self):
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index not in self.get_adsList()]]
        return tmpAtoms

    def get_optAtoms(self):
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index not in self.get_adsList() and\
               a.index not in self.get_bufList()]]
        return tmpAtoms

    def get_covalRadii(self):
        return [covalRadii[i] for i in self.get_allAtoms().numbers]

    def get_contactMat(self, scale=1.0):
        '''
        returns the standard covalent bondlength between all atoms
        '''
        cvRad = self.get_covalRadii()
        contact = np.tile(cvRad, [len(cvRad), 1])
        contact += contact.T
        contact *= scale
        np.fill_diagonal(contact, 0)
        return contact

    def get_allDistances(self):
        return self.get_allAtoms().get_all_distances(mic=True)

    def has_badContact(self, tolerance=0):
        diff = self.get_allDistances() - self.get_contactMat(scale=1-tolerance)
        return diff.min() < 0



    def set_allAtoms(self, newAllAtoms):
        newAllAtoms.wrap()
        self.allAtoms = newAllAtoms
        self.update()

    def set_allPos(self, newAllPos):
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.set_positions(newAllPos)
        self.set_allAtoms(tmpAtoms)

    def set_adsAtoms(self, newAdsAtoms):
        tmpAtoms = self.get_subAtoms()
        tmpAtoms.extend(newAdsAtoms)
        self.set_allAtoms(tmpAtoms)

    # def set_positions(self, newPos):
    #     tmpAtoms = self.get_allAtoms()
    #     tmpAtoms.set_positions(newPos)
    #     self.set_allAtoms(tmpAtoms)

    def get_chemical_symbols(self):
        return self.get_allAtoms().get_chemical_symbols().copy()

    def get_constraints(self):
        return self.subAtoms.constraints

    def get_topLayerList(self, depth = 1):
        allZ = self.get_pos[:,2]
        return [a.index for a in self.allAtoms\
            if max(allZ) - allZ[a.index] <= depth]

    def wrap(self):
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.wrap()
        self.set_allAtoms(tmpAtoms)

    def merge_adsorbate(self, adsAtoms):
        tmpAtoms = self.get_allAtoms()
        tmpAtoms.extend(adsAtoms)
        self.set_allAtoms(tmpAtoms)

    def remove_adatom(self, rmList=None):
        if rmList is None:
            rmList = list(range(len(oldAdsList)))
        if type(rmList) is not list:
            rmList = [rmList]
        adsList = [self.get_adsList()[i] for i in rmList]
        tmpAtoms = self.get_allAtoms()  
        del tmpAtoms[adsList]
        self.set_allAtoms(tmpAtoms)

    def rattle(self, stdev = 0.1, zEnhance=False):
        '''
        enhances the atoms with higher position
        '''
        tmpAtoms = self.get_allAtoms()
        pos = tmpAtoms.get_positions()
        rattleVec = np.random.normal(scale=stdev, size=pos.shape)
        if zEnhance:
            zBuf = self.get_bufAtoms().get_positions()[:,2]
            rattleVec = (rattleVec.T * (pos[:,2]-zBuf.min())/(pos[:,2].max()-zBuf.min())).T
        self.set_allPos(pos + rattleVec)

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
            trajectory=None,
            logfile=None
        )
        geomOpt.run(fmax = 0.01, steps = nsteps)
        print(' - L-J pre-optimization: RMSD = %.3f Angstroms'%geom.RMSD(self.get_allAtoms(), tmpAtoms))
        self.set_allAtoms(tmpAtoms)

    def preopt_hooke(self, cutoff = 1.5,
        toler=0.2, stepsize=0.05, nsteps=200):
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
        print(' - Hookean pre-optimization: RMSD = %.3f Angstroms'%geom.RMSD(self.get_allAtoms(), tmpAtoms))
        self.set_allAtoms(tmpAtoms)

    def write(self, fileName):
        fio.write(fileName, self.get_allAtoms())

    def draw(self, key='CPK', title=''):
        from gocia.utils import visualize
        if key == 'CPK':
            visualize.draw_CPKsurf(self, outName=self.tags, title=title)
        if key == 'BS':
            visualize.draw_BSsurf(self, outName=self.tags, title=title, pseudoBond=True)







        


