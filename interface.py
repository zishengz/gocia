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
from gocia import frag
from gocia.geom.build import grow_frag

import ase.io as fio
from ase.atoms import Atoms
from ase.constraints import FixAtoms
from ase.build.tools import sort
import numpy as np
import json

from re import A, I
import random
from scipy.spatial.transform import Rotation as R

class Interface:
    def __init__(self,
                 allAtoms=None,
                 subAtoms=None,
                 fixList=None,
                 bufList=None,
                 adsList=None,
                 fragList=None,
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
                self.allAtoms.wrap()
            else:
                self.allAtoms   = allAtoms

        self.fixList    = self.subAtoms.constraints[0].get_indices()
        self.bufList = [i for i in list(range(len(self.subAtoms)))\
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

        if fragList is not None:
            self.fragList = fragList
            self.set_fragList(self.fragList)
        else:
            self.fragList = None

    def __len__(self):
        return len(self.allAtoms)

    def copy(self):
        import copy
        myCopy = self.__class__(
            tags=copy.deepcopy(self.tags),
            subAtoms=copy.deepcopy(self.subAtoms),
            allAtoms=copy.deepcopy(self.allAtoms),
            zLim=copy.deepcopy(self.zLim),
            info=copy.deepcopy(self.info)
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
        self.tags = self.get_formula()
        self.info = self.allAtoms.info # Updates Interface info to match allAtoms info

    def sort(self):
        self.sortAds()

    def sortAds(self):
        ads = self.get_adsAtoms().copy()
        ads = sort(ads)
        self.set_adsAtoms(ads)

    def sortAds_frag(self, atoms=None):
        # Obtain atoms and list of fragments in adsorbate
        if atoms is not None:
            adsAtoms = atoms.copy()
            fragList = adsAtoms.info['adsorbate_fragments']
        else:
            adsAtoms = self.get_adsAtoms().copy()
            fragList = self.get_fragList()
        #print([self.get_allAtoms()[f].get_chemical_formula() for f in fragList])

        if not fragList:
            adsSort = adsAtoms
        else:
            # Make lists for atom symbols and fragment indices
            tags = adsAtoms.get_chemical_symbols()
            fragFlat = [atmId for frag in fragList for atmId in frag]
            # Sort fragment index list in same manner as sort atoms by chemical symbol
            ## Modified sort() provided by ase in tools.py
            deco = sorted([(tag, i) for i, tag in enumerate(tags)])
            indices = [i for tag, i in deco]
            adsSort = adsAtoms[indices]
            fragFlat = sorted(fragFlat) # JUST TRYING SOMETHING WILD THIS WORKED ! WHY NECESSARY? Which spot needed this? growMut I think
            fragSort = [fragFlat[i] for i in indices] #Problematic if partially sorted previouly? fragFlat not same as adsList!
            # Remake fragment list with new indices after sorting
            adsSort.info['adsorbate_fragments'] = frag.remake(fragList,fragSort,fragFlat)
        # either return sorted atoms object or set adsorbate atoms
        if atoms is not None:
            return adsSort
        else:
            self.set_adsAtoms_frag(adsSort) 
            return

    def print(self):
        print('#TAG: ', self.tags)
        print(' |-Fixed atoms:    ', self.get_fixAtoms().symbols)
        print(' |-Buffer atoms:   ', self.get_bufAtoms().symbols)
        if self.get_adsList() == []:
            print(' |-Adsorbates:     None')
        else:    
            print(' |-Adsorbates:     ',self.get_adsAtoms().symbols)
        if self.get_fragList() == []:
            print(' |-Fragments:     None')
        else:    
            print(' |-Fragment:     ',self.get_adsAtoms().symbols)
        print(' |-Buffer region:   Z = %.3f to %.3f'%\
            (self.zLim[0], self.zLim[1]))
        print(' |-Info:           ', self.info)
        if self.fragList is not None:
            print([self.get_allAtoms()[f].get_chemical_formula()for f in self.fragList])
        print('\n')

    def has_broken_frag(self):
        my_fragList = self.get_fragList()
        my_fragAtoms = [self.get_allAtoms()[f] for f in my_fragList]
        # Check connectivity
        flags = [len(geom.get_fragments(ads))==1 for ads in my_fragAtoms]
        if False in flags:
            return True
        else:
            return False

        
    def get_atomic_numbers(self):
        return self.allAtoms.get_atomic_numbers()

    def get_formula(self):
        return self.allAtoms.get_chemical_formula()

    def get_fragNames(self):
        fragNames = [self.get_allAtoms()[frag].get_chemical_formula() for frag in self.get_fragList()]
        return fragNames

    def get_cell(self):
        return self.cellParam.copy()

    def get_pbc(self):
        return self.pbcParam.copy()

    def get_pos(self):
        return self.get_allAtoms().get_positions().copy()

    def get_fixBufPos(self):
        return self.get_pos()[[i for i in list(range(len(self)))\
                    if i not in self.get_adsList()]]

    def get_subPos(self):
        return self.get_subAtoms().get_positions().copy()

    def get_fixList(self):
        return self.fixList.copy()

    def get_bufList(self):
        return self.bufList.copy()

    def get_adsList(self):
        return self.adsList.copy()

    def get_bridPos(self):
        return self.get_bridAtoms().get_positions().copy()

    def get_fragList(self):
        # gets list of fragments from interface info as updated when set_allAtoms
        if 'adsorbate_fragments' in self.info:
            return self.info['adsorbate_fragments']
        else:
            return []

    # Bridle atom is atom that steers fragment behavior
    def get_bridList(self): # assumes for now that the first atom is bridle atom
        return [frag[0] for frag in self.get_fragList()]

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

    def get_bridAtoms(self):
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index not in self.get_bridList()]]
        return tmpAtoms

    def get_optAtoms(self):
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index not in self.get_adsList() and\
               a.index not in self.get_bufList()]]
        return tmpAtoms

    def get_fixBufAtoms(self):
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index in self.get_adsList()]]
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

    def has_outsideBox(self):
        tmp = []
        ads = self.get_adsAtoms()
        for a in ads:
            if a.position[2] < min(self.zLim) or a.position[2] > max(self.zLim):
                tmp.append(a.index)
        if len(tmp) > 0:
            return True
        else:
            return False

    def get_outsideBox(self):
        # indeces are in adsAtoms list, not in allAtoms list!
        tmp = []
        ads = self.get_adsAtoms()
        for a in ads:
            if a.position[2] < min(self.zLim) or a.position[2] > max(self.zLim):
                tmp.append(a.index + len(self.get_allAtoms()) - len(self.get_adsAtoms()))
        return tmp

    def del_outsideBox(self):
        list_del = self.get_outsideBox()
        if len(list_del) > 0:
            print('Delete ', list_del)
            all = self.get_allAtoms()
            del all[list_del]
            self.set_allAtoms(all)

    def del_outsideBox_frag(self):
        list_del = self.get_outsideBox()
        list_frag = self.get_fragList()
        print(list_del)
        list_del_frag = []
        list_frag_new = []
        for f in list_frag:
            for d in list_del:
                if d in f:
                    list_del_frag += f
                    continue
            list_frag_new.append(
                [ii - len(list_del_frag) for ii in f]
                )
        print('list including whole fragments', list_del_frag)
        if len(list_del_frag) > 0:
            print('Delete ', list_del_frag)
            all = self.get_allAtoms()
            del all[list_del_frag]
            self.set_allAtoms(all)
            self.set_fragList(list_frag_new)

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

    def set_fixBufPos(self, newPos):
        tmpAtoms = self.get_subAtoms()
        tmpAtoms.set_positions(newPos)
        tmpAtoms.extend(self.get_adsAtoms())
        tmpAtoms.info = self.get_adsAtoms().info
        self.set_allAtoms(tmpAtoms)

    def set_adsAtoms(self, newAdsAtoms):
        tmpAtoms = self.get_fixBufAtoms()
        tmpAdsAtoms = newAdsAtoms.copy()
        if len(newAdsAtoms) > 0:
            tmpAdsAtoms = sort(newAdsAtoms)
        tmpAtoms.extend(tmpAdsAtoms)
        self.set_allAtoms(tmpAtoms)

    def set_adsAtoms_frag(self, newAdsAtoms):
        tmpAtoms = self.get_fixBufAtoms()
        tmpAdsAtoms = newAdsAtoms.copy()
        if len(newAdsAtoms) > 0:
            tmpAdsAtoms = self.sortAds_frag(newAdsAtoms) # sortAds tracks changes to fragList
        tmpAtoms.extend(tmpAdsAtoms)
        tmpAtoms.info = tmpAdsAtoms.info # all atoms take on info from sorted new adsorbate atoms
        self.set_allAtoms(tmpAtoms)
        self.set_fragList(tmpAdsAtoms.info['adsorbate_fragments'])

    def set_fragList(self, my_fragList):
        # sets the list of fragments as provided
        if type(my_fragList) is list:
            self.info['adsorbate_fragments'] = my_fragList
            self.fragList = my_fragList
        elif type(my_fragList) is str:
            try:
                self.info['adsorbate_fragments'] = eval(my_fragList)
                self.fragList = eval(my_fragList)
            except:
                print('fragLIst must be list or str!')

    # def set_positions(self, newPos):
    #     tmpAtoms = self.get_allAtoms()
    #     tmpAtoms.set_positions(newPos)
    #     self.set_allAtoms(tmpAtoms)

    def get_chemical_symbols(self):
        return self.get_allAtoms().get_chemical_symbols().copy()

    def get_constraints(self):
        return self.subAtoms.constraints

    def get_topLayerList(self, depth = 1):
        allZ = self.get_pos()[:,2]
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

    def remove_adatom_old(self, rmList=None):
        if rmList is None:
            rmList = list(range(len(self.get_adsList())))
        if type(rmList) is not list:
            rmList = [rmList]
        adsList = [self.get_adsList()[i] for i in rmList]
        tmpAtoms = self.get_allAtoms()  
        del tmpAtoms[adsList]
        self.set_allAtoms(tmpAtoms)

    def add_adsFrag(self, adsAtoms): # assumes adding single adsorbate molecule
        tmpAtoms = self.get_allAtoms().copy()
        oldAdsList = self.get_adsList().copy()
        tmpAtoms.extend(adsAtoms)
        # length of all atoms is the sum of substrate, previously added adsorbate, and new adsorbate atoms
        # (calculating it this way enables not having to set_allAtoms() then get_adsList())
        allLen = len(self.get_subAtoms()) + len(self.get_adsList()) + len(adsAtoms)
        newAdsList = [i for i in list(range(allLen))\
            if i not in list(range(len(self.get_subAtoms())))]  
        fragList = []
        if 'adsorbate_fragments' in tmpAtoms.info:
            fragList = tmpAtoms.info['adsorbate_fragments']
        diff = list(set(newAdsList) - set(oldAdsList))
        # only append if not already in list as otherwise keep adding same fragment while badStructure
        if diff not in fragList:
            fragList.append(diff)
        tmpAtoms.info['adsorbate_fragments'] = fragList
        self.set_allAtoms(tmpAtoms)

    def remove_adsFrag(self, rmList=None):       #currently only able to remove one adsorbate at a time
        fragList = self.get_fragList()
        rmFragList = []
        if not fragList:
            print("WARNING: No adsorbate fragment to remove.")
            return
        # properly choose and format fragments to remove as rmList can vary
        if rmList is None: 
            rmFragList = [fragList[-1]] # if not rmList provided, just remove the last fragment
        if type(rmList) is list:
            if all(isinstance(x,list) for x in rmList) and all(isinstance(y,int) for x in rmList for y in x):
                rmFragList = rmList # if rmList is 
            if all(isinstance(x,str) for x in rmList):
                print('WARNING: Currently unable to remove fragments inputed as list of strings.')
            if all(isinstance(x,int) for x in rmList):
                rmFragList = [rmList]
        else:
            if type(rmList) is str:
                rand = random.randint(0,len(fragList)-1)
                rmFragList = [fragList[rand]]
        # Got the fragment in a list so now need to remove it from fragList and update adsAtoms
        rmFlat = [atmId for frag in rmFragList for atmId in frag]
        tmpAtoms = self.get_allAtoms()
        del tmpAtoms[[a.index for a in tmpAtoms\
            if a.index in rmFlat]]
        delFragList = [frag for frag in fragList if frag not in rmFragList]
        if not delFragList:
            tmpAtoms.info['adsorbate_fragments'] = []
        else:
            fragFlatSort = frag.flatsort(delFragList)
            reindexList = list(range(len(self.get_subAtoms()),int(len(self.get_subAtoms())+len(fragFlatSort))))
            tmpAtoms.info['adsorbate_fragments'] = frag.remake(delFragList,fragFlatSort,reindexList)
        self.set_allAtoms(tmpAtoms) # Do we need to sort before setting? Currently it keeps the same order just with fragments removed

    def detect_fragList(self, scale=1.0, update = False):
        ads = self.get_adsAtoms()
        n_subs = len(self.get_subAtoms())
        my_fragList = geom.get_fragments(ads, scale=scale)
        my_fragList = [[i+n_subs for i in f] for f in my_fragList]
        if update:
            self.set_fragList(my_fragList)
        return my_fragList

    def rattle(self, stdev = 0.1, zEnhance=False):
        '''
        enhances the atoms with higher position
        '''
        tmpAtoms = self.get_allAtoms()
        pos = tmpAtoms.get_positions()
        zBuf = self.get_bufAtoms().get_positions()[:,2]
        rattleVec = np.random.normal(scale=stdev, size=pos.shape)
        if zEnhance and pos[:,2].max()-zBuf.min() != 0:
            rattleVec = (rattleVec.T * (pos[:,2]-zBuf.min())/(pos[:,2].max()-zBuf.min())).T
        self.set_allPos(pos + rattleVec)

    def rattleMut(self, stdev = 0.25, mutRate = 0.5, zEnhance=True):
        '''
        enhances the atoms with higher position
        '''
        print(' |- Rattle mutation!')
        tmpAtoms = self.get_allAtoms()
        pos = tmpAtoms.get_positions()
        zBuf = self.get_bufAtoms().get_positions()[:,2]
        rattleVec = np.random.normal(scale=stdev, size=pos.shape)
        if zEnhance and pos[:,2].max()-zBuf.min() != 0:
            rattleVec = (rattleVec.T * (pos[:,2]-zBuf.min())/(pos[:,2].max()-zBuf.min())).T
        for i in self.get_bufList():
            if np.random.rand() < mutRate:
                pos[i] += rattleVec[i] * 0.5
                if pos[i][2] > max(zBuf): pos[i][2] = max(zBuf)
                if pos[i][2] < min(zBuf): pos[i][2] = min(zBuf)
        for i in self.get_adsList():
            if np.random.rand() < mutRate:
                pos[i] += rattleVec[i]
                if pos[i][2] > max(self.zLim): pos[i][2] = max(self.zLim)
                if pos[i][2] < min(self.zLim): pos[i][2] = min(self.zLim)
        self.set_allPos(pos)
    
    def rattleMut_buffer(self, stdev = 0.25, mutRate = 0.5, zEnhance=True):
        print(' |- Rattle mutation! -- buffer atoms')
        tmpAtoms = self.get_allAtoms()
        pos = tmpAtoms.get_positions()
        zBuf = self.get_bufAtoms().get_positions()[:,2]
        rattleVec = np.random.normal(scale=stdev, size=pos.shape)
        if zEnhance and pos[:,2].max()-zBuf.min() != 0:
            rattleVec = (rattleVec.T * (pos[:,2]-zBuf.min())/(pos[:,2].max()-zBuf.min())).T
        for i in self.get_bufList():
            if np.random.rand() < mutRate:
                pos[i] += rattleVec[i]
                if pos[i][2] > max(zBuf): pos[i][2] = max(zBuf)
                if pos[i][2] < min(zBuf): pos[i][2] = min(zBuf)
        self.set_allPos(pos)


    def rattleMut_frag(self, stdev = 0.2, mutRate = 0.5, zEnhance=False, toler=0.5):

        # Initialize as before
        print(' |- Rattle mutation! -- fragments')
        tmpAtoms = self.get_allAtoms()
        pos = tmpAtoms.get_positions()
        zBuf = self.get_bufAtoms().get_positions()[:,2]
        rattleVec = np.random.normal(scale=stdev, size=pos.shape)
        # Enhance the rattling of atoms with higher position
        if zEnhance and pos[:,2].max()-zBuf.min() != 0:
            rattleVec = (rattleVec.T * (pos[:,2]-zBuf.min())/(pos[:,2].max()-zBuf.min())).T
        # Rattle buffer atoms as before, shifting position by rattleVec and ensuring within z-limits
        for i in self.get_bufList():
            if np.random.rand() < mutRate:
                pos[i] += rattleVec[i] * 0.5
                if pos[i][2] > max(zBuf): pos[i][2] = max(zBuf)
                if pos[i][2] < min(zBuf): pos[i][2] = min(zBuf)

        # Now for rattling fragments . . . 
        fragList = self.get_fragList()
        fragNames = self.get_fragNames()
        bridList = self.get_bridList()
        for fragID in range(len(fragList)):
            # Rattle fragments that are single atoms
            if len(fragList[fragID]) == 1 and np.random.rand() < mutRate:
                i = fragList[fragID][0]
                pos[i] += rattleVec[i]
                if pos[i][2] > max(self.zLim): pos[i][2] = max(self.zLim)
                if pos[i][2] < min(self.zLim): pos[i][2] = min(self.zLim) 
            # Rotate multiatom fragments as a whole unit
            if len(fragList[fragID]) > 1 and np.random.rand() < mutRate:
                # Need name of fragment
                fragName = fragNames[fragID]
                # Rotate fragment positions until no collisions
                keepStructure = False
                while not keepStructure: 
                    # Position of bridle atom is reference for positioning whole fragment
                    bridPos = pos[bridList[fragID]]
                    # Now we rattle bridle atom position
                    bridPos += rattleVec[bridList[fragID]]
                    # Make nicely oriented atoms object
                    fragTemp = ''
                    if fragName == 'CO':
                        fragTemp = Atoms('CO',[(0, 0, 0),(0, 0, 1.43)])
                    elif fragName == 'H':
                        fragTemp = Atoms('H',[(0,0,0)])
                    elif fragName == 'H2O':
                        fragTemp = Atoms('OH2',[(0,0,0),(0.758602,0,0.504284),(-0.758602,0,0.504284)])
                    else:
                        print('Unknown fragment: must add option to grow {} in twistMut() in interface.py'.format(fragName))
                    # Ensure fragTemp sorted in same way as order of fragment atoms in pos
                    fragTemp = sort(fragTemp)
                    # Get well-ordered positions
                    posTemp = fragTemp.get_positions()
                    # Rotate positions about z-axis
                    rot = R.from_rotvec(- random.random()*np.pi * np.array([0, 0, 1]))
                    posTemp = rot.apply(posTemp)
                    # Rotate positions by random angle (less than 45 deg) about random axis 
                    rot = R.from_rotvec(- random.random()*np.pi/4 * geom.rand_direction())
                    posTemp = rot.apply(posTemp)
                    # Set positions relative to rattled bridle atom position
                    posNew = posTemp + bridPos
                    # Update positions
                    for i in range(len(fragList[fragID])):
                        pos[fragList[fragID][i]] = posNew[i]
                    # Keep structure unless it has bad contacts (might not even be necessary as didn't check before?)
                    tmpTest = self.copy()
                    tmpTest.set_allPos(pos)
                    if not tmpTest.has_badContact(tolerance=toler): # Hmm do I need to check for bad contact in the buffer atoms?
                        keepStructure = True 
        # Set positions
        self.set_allPos(pos)

    def transMut(self, transVec=[[-2,2],[-2,2]]):
        tmpAds = self.get_adsAtoms()
        myAxis = np.random.choice([0,1],size=1)[0]
        #example of the translational vector: [[-3, 3],[-2, 2]]
        vecPeriod = np.random.choice(transVec[myAxis],size=1)[0]
        print(' |- Translation mutation! axis=%i, periodicity=%i'%(myAxis, vecPeriod))
        tmpAds.set_positions(tmpAds.get_positions()+tmpAds.get_cell()[myAxis]/vecPeriod)
        tmpAds.wrap()
        self.set_adsAtoms(tmpAds)

    def permuteMut(self):
        tmpAds = self.get_adsAtoms()
        myCell = np.array(tmpAds.get_cell())
        adsCom = tmpAds.get_center_of_mass()
        adsCom = geom.cart2frac(adsCom, myCell)
#        adsCom = np.random.rand(3)
        myAxis = np.random.choice([0,1],size=1)[0]
        myBase = np.random.choice([0,1],size=1)[0]
        print(' |- Permutation mutation!')
        tmpAdsHalf = tmpAds.copy()
        tmpAdsFrac = geom.cart2frac(tmpAdsHalf.get_positions(), myCell)
        if myBase == 0:
            del tmpAds[[i for i in range(len(tmpAdsHalf))\
                if tmpAdsFrac[i][myAxis] > adsCom[myAxis] ]]
        else:
            del tmpAds[[i for i in range(len(tmpAdsHalf))\
                if tmpAdsFrac[i][myAxis] < adsCom[myAxis] ]]
        tmpTmp = tmpAds.copy()
        tmpTmp.set_positions(tmpTmp.get_positions() + myCell[myAxis]/2)
        tmpAds.extend(tmpTmp)
        tmpAds.wrap()
        self.set_adsAtoms(tmpAds)

    def permuteMut_frag_old(self):
        # Obtain adsorbate atoms and their center of mass within the unit cell
        tmpAds = self.get_adsAtoms()
        myCell = np.array(tmpAds.get_cell())
        adsCom = tmpAds.get_center_of_mass()
        adsCom = geom.cart2frac(adsCom, myCell)
        # Choose parameters designating two quadrants to permute
        myAxis = np.random.choice([0,1],size=1)[0]
        myBase = np.random.choice([0,1],size=1)[0]
        print(' |- Permutation mutation!')
        # Reindex frags from All to Ads: [[144,146],[145,147]] to [[0, 2], [1, 3]] 
        fragList_tmpAds = frag.remake(self.get_fragList(),self.get_adsList(),[a - len(self.get_subAtoms()) for a in self.get_adsList()])
        # Get atoms which bridle fragments and their relative positions
        bridleAtms = [a - len(self.get_subAtoms()) for a in self.get_bridList()] # length of bridleList should equal length of fragList_tmpAds
        tmpAdsFrac = geom.cart2frac(tmpAds.get_positions(), myCell)
        # Remove fragments which aren't bridled within choosen quadrants
        if myBase == 0:
            del tmpAds[[i for f in range(len(fragList_tmpAds)) for i in fragList_tmpAds[f] if tmpAdsFrac[bridleAtms[f]][myAxis] > adsCom[myAxis] ]]
            fragList_tmpAds = [fragList_tmpAds[f] for f in range(len(fragList_tmpAds)) if not tmpAdsFrac[bridleAtms[f]][myAxis] > adsCom[myAxis] ]
        else: 
            del tmpAds[[i for f in range(len(fragList_tmpAds)) for i in fragList_tmpAds[f] if tmpAdsFrac[bridleAtms[f]][myAxis] < adsCom[myAxis] ]]
            fragList_tmpAds = [fragList_tmpAds[f] for f in range(len(fragList_tmpAds)) if not tmpAdsFrac[bridleAtms[f]][myAxis] < adsCom[myAxis] ]
        # Reindex remaining fragments 
        if not len(fragList_tmpAds) > 0:
            fragList_tmpAds = []
        else:
            fragList_tmpAds = frag.remake(fragList_tmpAds,frag.flatsort(fragList_tmpAds),range(len(tmpAds)))
        # Duplicate and translate adsorbate atoms
        tmpTmp = tmpAds.copy()
        tmpTmp.set_positions(tmpTmp.get_positions() + myCell[myAxis]/2)
        tmpAds.extend(tmpTmp) 
        # Make fragment list including originals and duplicates 
        flat = frag.flatten(fragList_tmpAds)
        fill = flat + [len(flat) + i for frag in fragList_tmpAds for i in frag]
        fragList_dup = fragList_tmpAds*2
        tmpAds.info['adsorbate_fragments'] = frag.refill(fragList_dup,fill)
        # Clean-up adsorbate atoms
        tmpAds = self.sortAds_frag(tmpAds)
        tmpAds.wrap()
        # Reindex frags from Ads to All: [[0, 2], [1, 3]] to [[144,146],[145,147]] 
        tmpAdsInd = [i for i in list(range(len(tmpAds)))]
        tmpAds.info['adsorbate_fragments'] = frag.remake(tmpAds.info['adsorbate_fragments'], tmpAdsInd, [a + len(self.get_subAtoms()) for a in tmpAdsInd])
        self.set_adsAtoms(tmpAds)

    def permuteMut_frag(self):
        # Obtain adsorbate atoms and their center of mass within the unit cell
        tmpAds = self.get_adsAtoms()
        myCell = np.array(tmpAds.get_cell())
        adsCom = tmpAds.get_center_of_mass()
        adsCom = geom.cart2frac(adsCom, myCell)
        # Choose parameters designating two quadrants to permute
        myAxis = np.random.choice([0,1],size=1)[0]
        myBase = np.random.choice([0,1],size=1)[0]
        print(' |- Permutation mutation!')
        # Reindex frags from All to Ads: [[144,146],[145,147]] to [[0, 2], [1, 3]] 
        fragList_tmpAds = frag.remake(self.get_fragList(),self.get_adsList(),[a - len(self.get_subAtoms()) for a in self.get_adsList()])
        # Get atoms which bridle fragments and their relative positions
        bridleAtms = [a - len(self.get_subAtoms()) for a in self.get_bridList()] # length of bridleList should equal length of fragList_tmpAds
        tmpAdsFrac = geom.cart2frac(tmpAds.get_positions(), myCell)
        # Remove fragments which aren't bridled within choosen quadrants
        if myBase == 0:
            del tmpAds[[i for f in range(len(fragList_tmpAds)) for i in fragList_tmpAds[f] if tmpAdsFrac[bridleAtms[f]][myAxis] > adsCom[myAxis] ]]
            fragList_tmpAds = [fragList_tmpAds[f] for f in range(len(fragList_tmpAds)) if not tmpAdsFrac[bridleAtms[f]][myAxis] > adsCom[myAxis] ]
        else: 
            del tmpAds[[i for f in range(len(fragList_tmpAds)) for i in fragList_tmpAds[f] if tmpAdsFrac[bridleAtms[f]][myAxis] < adsCom[myAxis] ]]
            fragList_tmpAds = [fragList_tmpAds[f] for f in range(len(fragList_tmpAds)) if not tmpAdsFrac[bridleAtms[f]][myAxis] < adsCom[myAxis] ]
        # Reindex remaining fragments 
        if not len(fragList_tmpAds) > 0:
            fragList_tmpAds = []
        else:
            fragList_tmpAds = frag.remake(fragList_tmpAds,frag.flatsort(fragList_tmpAds),range(len(tmpAds)))
        # Duplicate and translate adsorbate atoms
        tmpTmp = tmpAds.copy()
        tmpTmp.set_positions(tmpTmp.get_positions() + myCell[myAxis]/2)
        tmpAds.extend(tmpTmp) 
        # Make fragment list including originals and duplicates 
        flat = frag.flatten(fragList_tmpAds)
        fill = flat + [len(flat) + i for frag in fragList_tmpAds for i in frag]
        fragList_dup = fragList_tmpAds*2
        fragList_dup = frag.refill(fragList_dup,fill)
        #fragList_dup = frag.transposeUp(fragList_dup, self)
        tmpAds.info['adsorbate_fragments'] = fragList_dup
        self.set_fragList(fragList_dup)
        print(fragList_dup)
        # Clean-up adsorbate atoms
        tmpAds = self.sortAds_frag(tmpAds)
        tmpAds.wrap()
        # Reindex frags from Ads to All: [[0, 2], [1, 3]] to [[144,146],[145,147]] 
        tmpAdsInd = [i for i in list(range(len(tmpAds)))]
        fragList_final = frag.remake(tmpAds.info['adsorbate_fragments'], tmpAdsInd, [a + len(self.get_subAtoms()) for a in tmpAdsInd])
        tmpAds.info['adsorbate_fragments'] = fragList_final
        self.set_adsAtoms(tmpAds)
        self.set_fragList(fragList_final)

    def leachMut(self, elemList):
        print(' |- Leaching mutation:', end = '\t')
        tmpAds = self.get_adsAtoms()
        nads = len(tmpAds)
        while len(tmpAds) == nads:
            myDel = np.random.choice(list(range(nads)),size=1)[0]
            if tmpAds.get_chemical_symbols()[myDel] in elemList:
                print(tmpAds.get_chemical_symbols()[myDel])
                del tmpAds[myDel]
        self.set_adsAtoms(tmpAds)

    def leachMut_frag(self, fragPool):
        print(' |- Leaching mutation:', end = '\t')
        fragList = self.get_fragList()
        nfrags = len(fragList)
        while len(self.get_fragList()) == nfrags:
            myDel = np.random.choice(list(range(nfrags)),size=1)[0]
            print(self.get_fragNames()[myDel])
            if self.get_fragNames()[myDel] in fragPool:
                print(self.get_fragNames()[myDel])
                self.remove_adsFrag(fragList[myDel])

    def growMut(self, elemList):
        print(' |- Growth mutation:', end = '\t')
        from gocia.geom.build import grow_adatom
        tmpInterfc = self.copy()
        myElem = np.random.choice(elemList, size=1)[0]
        print(myElem)
        tmpInterfc = grow_adatom(
            tmpInterfc,
            myElem,
            zLim = self.zLim,
        )
        self.set_allAtoms(tmpInterfc.get_allAtoms())

    def growMut_frag(self, fragPool):
        print(' |- Growth mutation:', end = '\t')
        from gocia.geom.build import grow_frag
        tmpInterfc = self.copy()
        myFrag = np.random.choice(fragPool, size=1)[0]
        tmpInterfc = grow_frag(
            tmpInterfc,
            [myFrag],
            zLim = self.zLim,
        )
        self.set_allAtoms(tmpInterfc.get_allAtoms())

    def moveMut_frag(self, fragPool):
        # DO NOT USE IT FOR NOW!
        # The popGrandCanonPoly has a better implementation which includes constraints
        print(' |- Move mutation (Leach & Grow)', end = '\n')
        tmpInterfc = self.copy()
        myFrag = np.random.choice(fragPool, size=1)[0]
        tmpInterfc.leachMut_frag([myFrag])
        tmpInterfc.growMut_frag([myFrag])
        self.set_allAtoms(tmpInterfc.get_allAtoms())

    def growMut_box(self, elemList, xyzLims, bondRejList = None, constrainTop=False):
        print(' |- Growth mutation:', end = '\t')
        from gocia.geom.build import boxSample_adatom
        tmpInterfc = self.copy()
        myElem = np.random.choice(elemList, size=1)[0]
        print(myElem)
        tmpInterfc = boxSample_adatom(
            tmpInterfc,
            myElem,
            xyzLims=xyzLims,
            bondRejList=bondRejList,
            constrainTop=constrainTop
        )
        if tmpInterfc is not None:
            self.set_allAtoms(tmpInterfc.get_allAtoms())

    def growMut_box_frag(self, fragPool, xyzLims, bondRejList = None, constrainTop=False):
        print(' |- Growth mutation:', end = '\t')
        tmpInterfc = self.copy()
        myFrag = np.random.choice(fragPool, size=1)[0]
        print(myFrag)
        from gocia.geom.build import boxSample_frag
        tmpInterfc = boxSample_frag(
            tmpInterfc,
            [myFrag],
            xyzLims=xyzLims,
            bondRejList=bondRejList,
            constrainTop=constrainTop
        )
        if tmpInterfc is not None:
            self.set_allAtoms(tmpInterfc.get_allAtoms())

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
            visualize.draw_BSsurf(self, outName=self.tags, title=title, hBond=True)
            #visualize.draw_BSsurf(self, outName=self.tags, title=title, pseudoBond=True)







        


