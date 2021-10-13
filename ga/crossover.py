import sys, os
import numpy as np
import ase.io as fio
from ase.atoms import Atoms
from gocia.interface import Interface

def dist2lin(p1, p2, p3):
    '''
    p1 and p2 defines the line
    returns the distance with signs (+/-)
    '''
    return np.cross(p2-p1, p1-p3)/np.linalg.norm(p2-p1)

def split_2d(atoms1, atoms2):
    '''
    Returns the center of the COMs of the parents
        and a random slope in x-y plane.
    '''
    com1 = atoms1.get_center_of_mass()
    com2 = atoms2.get_center_of_mass()
    direction = np.array([1,1])
    while sum(direction**2) >= 1:
        direction = 2*np.random.rand(2)-1
    com_share = (com1 + com2)[:2]/2
    return com_share, direction

def splice_2d(surf1, badPos, goodPos, center, direction):
    tmpSurf = surf1.copy()
    childPos = []
    for i in range(len(tmpSurf)):
        if len(goodPos[i]) == 1:
            childPos.append(goodPos[i][0])
        elif len(goodPos[i]) == 0:
            # Then there must be 2 in bad position list
            # We keep the one closer to splitline
            tmpDist = [dist2lin(center, center+direction, p[:2])\
                for p in badPos[i]]
            childPos.append(badPos[i][tmpDist.index(min(tmpDist))])
        elif len(goodPos[i]) == 2:
            # We keep the one closer to splitline
            tmpDist = [dist2lin(center, center+direction, p[:2])\
                for p in goodPos[i]]
            childPos.append(goodPos[i][tmpDist.index(min(tmpDist))])
    tmpSurf.set_allPos(childPos)
    return tmpSurf

def crossover_snsSurf_2d(surf1, surf2, tolerance=0.5):
    matPos, patPos = surf1.get_pos(), surf2.get_pos()
    goodPos, badPos = [], []
    isBADSTRUCTURE = True
    n_trial = 0
    while isBADSTRUCTURE:
        center, direction = split_2d(surf1.get_adsAtoms(), surf2.get_adsAtoms())
        matDist = [dist2lin(center, center+direction, i) for i in matPos[:, :2]]
        patDist = [dist2lin(center, center+direction, i) for i in patPos[:, :2]]
        n_trial += 1
        for i in range(len(surf1)):
            tmpGood, tmpBad = [], []
            if matDist[i] <= 0:
                tmpGood.append(matPos[i])
            else:
                tmpBad.append(matPos[i])
            if patDist[i] > 0:
                tmpGood.append(patPos[i])
            else:
                tmpBad.append(patPos[i])
            goodPos.append(tmpGood)
            badPos.append(tmpBad)
        # print(sum([len(l) for l in goodPos]))
        # print(sum([len(l) for l in badPos]))
        childSurf = splice_2d(surf1, badPos, goodPos, center, direction)
        isBADSTRUCTURE = childSurf.has_badContact(tolerance=tolerance)
        if n_trial > 100:
            isBADSTRUCTURE = False
            childSurf = None
    if childSurf is not None:
        print('\nOffspring is created at attempt #%i\t|Tolerance = %.3f'%\
            (n_trial, tolerance))
    return childSurf

def splice_2d_GC(surf1, badPos, goodPos, center, direction):
    tmpSurf = surf1.copy()
    childPos = []
    for i in range(len(tmpSurf)):
        if len(goodPos[i]) == 1:
            childPos.append(goodPos[i][0])
        elif len(goodPos[i]) == 0:
            # Then there must be 2 in bad position list
            # We keep the one closer to splitline
            tmpDist = [dist2lin(center, center+direction, p[:2])\
                       for p in badPos[i]]
            childPos.append(badPos[i][tmpDist.index(min(tmpDist))])
        elif len(goodPos[i]) == 2:
            # We keep the one closer to splitline
            tmpDist = [dist2lin(center, center+direction, p[:2])\
                       for p in goodPos[i]]
            childPos.append(goodPos[i][tmpDist.index(min(tmpDist))])
    tmpSurf.set_allPos(childPos)
    return tmpSurf

def crossover_snsSurf_2d_GC(surf1, surf2, tolerance=0.5):
    matPos, patPos = surf1.get_fixBufPos(), surf2.get_fixBufPos()
    matAds, patAds = surf1.get_adsAtoms(), surf2.get_adsAtoms()
    goodPos, badPos = [], []
    isBADSTRUCTURE = True
    n_trial = 0
    while isBADSTRUCTURE:
        n_trial += 1
        center, direction = split_2d(surf1.get_adsAtoms(), surf2.get_adsAtoms())

        # keep number of atoms constant in the buffer region!
        matDist = [dist2lin(center, center+direction, i) for i in matPos[:, :2]]
        patDist = [dist2lin(center, center+direction, i) for i in patPos[:, :2]]
        for i in range(len(matPos)):
            tmpGood, tmpBad = [], []
            if matDist[i] <= 0:
                tmpGood.append(matPos[i])
            else:
                tmpBad.append(matPos[i])
            if patDist[i] > 0:
                tmpGood.append(patPos[i])
            else:
                tmpBad.append(patPos[i])
            goodPos.append(tmpGood)
            badPos.append(tmpBad)
        childPos = []
        newSurf = surf1.copy()
        for i in range(len(matPos)):
            if len(goodPos[i]) == 1:
                childPos.append(goodPos[i][0])
            elif len(goodPos[i]) == 0:
                # Then there must be 2 in bad position list
                # We keep the one closer to splitline
                tmpDist = [dist2lin(center, center+direction, p[:2])\
                        for p in badPos[i]]
                childPos.append(badPos[i][tmpDist.index(min(tmpDist))])
            elif len(goodPos[i]) == 2:
                # We keep the one closer to splitline
                tmpDist = [dist2lin(center, center+direction, p[:2])\
                        for p in goodPos[i]]
                childPos.append(goodPos[i][tmpDist.index(min(tmpDist))])
        fixBufPos = childPos.copy()
        newSurf.set_fixBufPos(fixBufPos)

        # Then we work on the ads atoms
        matDist = [dist2lin(center, center+direction, i) for i in matAds.get_positions()[:, :2]]
        patDist = [dist2lin(center, center+direction, i) for i in patAds.get_positions()[:, :2]]
        newAdsElem, newAdsPos = [], []
#        print(matDist, patDist)
        for i in range(len(matDist)):
            if matDist[i] <= 0:
                newAdsElem.append(matAds.get_chemical_symbols()[i])
                newAdsPos.append(matAds.get_positions()[i])
        for i in range(len(patDist)):
            if patDist[i] > 0:
                newAdsElem.append(patAds.get_chemical_symbols()[i])
                newAdsPos.append(patAds.get_positions()[i])
        newSurf.set_adsAtoms(Atoms(newAdsElem, newAdsPos))
        newSurf.wrap()
        isBADSTRUCTURE = newSurf.has_badContact(tolerance=tolerance)
        if n_trial > 1000:
            isBADSTRUCTURE = False
            newSurf = None
    if newSurf is not None:
        print('\nOffspring is created at attempt #%i\t|Tolerance = %.3f'%\
            (n_trial, tolerance))
    return newSurf


        # print(sum([len(l) for l in goodPos]))
        # print(sum([len(l) for l in badPos]))
#        childSurf = splice_2d_GC(surf1, badPos, goodPos, center, direction)
    #     isBADSTRUCTURE = childSurf.has_badContact(tolerance=tolerance)
    #     if n_trial > 100:
    #         isBADSTRUCTURE = False
    #         childSurf = None
    # if childSurf is not None:
    #     print('\nOffspring is created at attempt #%i\t|Tolerance = %.3f'%\
    #         (n_trial, tolerance))
    # return childSurf