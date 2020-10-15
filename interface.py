#
# Created on Wed Oct 14 2020
#
# Copyright (c) 2020 Zisheng Zhang
# Alexandrova Lab, UCLA, United States

"""
This module defines the Surface object. 
"""

from gaia.data import elemSymbol, covalRadii
import numpy as np

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
            self.subAtoms   =subAtoms
            self.constraints= subAtoms.constraints
            self.fixList    = self.constraints[0].get_indices()
            self.bufferList = [i for i in list(range(len(subAtoms)))\
                if i not in self.fixList]
            self.cellParam  = allAtoms.get_cell()
            self.pbcParam   = allAtoms.get_pbc()            
        
        if allAtoms is not None:
            self.allAtoms   = allAtoms
            self.allPos     = self.allAtoms.get_positions()

        if zLim is None:
            self.zLim       = [min(self.allPos[:,2][self.bufferList]) + 0.1,
                    max(self.allPos[:,2][self.bufferList]) + 2]

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
        self.allAtoms = newAllAtoms

    def get_cell(self):
        return self.cellParam

    def get_adsList(self):
        return np.array(
                [[self.allAtoms.get_atomic_numbers()[i], i]\
                for i in list(range(len(self.allAtoms)))\
                if i not in list(range(len(self.subAtoms)))]
            )

    def get_topLayer(self, depth = 1):
        allZ = self.allPos[:,2]
        return [a.index for a in self.allAtoms\
            if max(allZ) - allZ[a.index] <= depth]


    def merge_adsorbate(self, adsAtoms):
        self.set_allAtoms(
            self.get_allAtoms().extend(adsAtoms)
        )
        self.update()

    def update(self):
        self.allPos  = self.allAtoms.get_positions()
        self.adsList = self.get_adsList()

    def print_info(self):
        print(self.tags)
        print(self.allAtoms)
        print('Buffer region:\tZ = %.3f to %.3f'%\
            (self.zLim[0], self.zLim[1]))
        print('Adsorbate:\t'+\
            str([elemSymbol[i] for i in self.adsList[:,0]]))




        


