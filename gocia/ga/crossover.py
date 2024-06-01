import sys, os
import numpy as np
import ase.io as fio
from ase.atoms import Atoms
from gocia.interface import Interface
from gocia import frag
from gocia import geom
from ase.build.tools import sort

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
    matPos, patPos = surf1.get_fixBufPos(), surf2.get_fixBufPos() # positions of fixed and buffer atoms in substrate
    matAds, patAds = surf1.get_adsAtoms(), surf2.get_adsAtoms() # adsorbate atoms
    goodPos, badPos = [], []
    # look until we find a viable structure
    isBADSTRUCTURE = True
    n_trial = 0
    while isBADSTRUCTURE:
        n_trial += 1
        # find center of mass and generate a random direction to define a line for splitting
        center, direction = split_2d(surf1.get_adsAtoms(), surf2.get_adsAtoms())

        # make list of distances between the position of every mother/father atom in the x,y-plane and the splitting line
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
        # make a child by taking proper atoms
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
        if n_trial > 100:
            isBADSTRUCTURE = False
            newSurf = None
    if newSurf is not None:
        print('\nOffspring is created at attempt #%i\t|Tolerance = %.3f'%\
            (n_trial, tolerance))
    return newSurf    

def crossover_snsSurf_2d_GC_poly(surf1, surf2, tolerance=0.5, bondRejList=None):
    # Get relevant positions from mother and father surfaces
    matFixBufPos, patFixBufPos = surf1.get_fixBufPos(), surf2.get_fixBufPos() # positions of fixed and buffer atoms, or equivalently, substrate atoms with get_subPos()
    matBridPos, patBridPos = surf1.get_bridPos(), surf2.get_bridPos() # positions of bridle atoms
    matAds, patAds = surf1.get_adsAtoms(), surf2.get_adsAtoms() # adsorbate atoms

    # Look until we find a viable child structure
    isBADSTRUCTURE = True
    n_trial = 0
    goodPos, badPos = [], []
    while isBADSTRUCTURE:
        n_trial += 1

        # Find center of mass and generate a random direction to define line for splitting
        center, direction = split_2d(matAds, patAds)

        ### Start with fixed and buffer atoms
        # Make list of distances between the position of every mother/father atom in the x,y-plane and the splitting line
        matDist = [dist2lin(center, center+direction, i) for i in matFixBufPos[:, :2]]
        patDist = [dist2lin(center, center+direction, i) for i in patFixBufPos[:, :2]]

        # Sort all mother/father atoms into good and bad
        # REQUIRES SAME NUMBER OF FIXED AND BUFFER ATOMS IN MOTHER AND FATHER
        for i in range(len(matFixBufPos)): 
            tmpGood, tmpBad = [], []
            if matDist[i] <= 0:     
                tmpGood.append(matFixBufPos[i])
            else:                  
                tmpBad.append(matFixBufPos[i])
            if patDist[i] > 0:    
                tmpGood.append(patFixBufPos[i])
            else:
                tmpBad.append(patFixBufPos[i])
            goodPos.append(tmpGood)
            badPos.append(tmpBad)

        # Make child by taking proper atoms
        childPos = []
        childSurf = surf1.copy()
        for i in range(len(matFixBufPos)):
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
        childSurf.set_fixBufPos(childPos)

        ### Then work on adsorbate fragment atoms
        # Take distance from bridle atom of each fragment to splitting line
        matDist = [dist2lin(center, center+direction, i) for i in matBridPos[:, :2]]
        patDist = [dist2lin(center, center+direction, i) for i in patBridPos[:, :2]]

        # !!! Similar strategy in permuteMut(): slide index, separately select and sort, fuse together, slide back
        # Get mother/father fragment lists and transpose indices down to only adsorbate atoms
        matFragList = frag.transposeDown(surf1.get_fragList(), surf1)
        patFragList = frag.transposeDown(surf2.get_fragList(), surf2)

        # print('matFragList: ', matFragList, [matAds[f].get_chemical_formula() for f in matFragList])
        # print('patFragList: ', patFragList, [patAds[f].get_chemical_formula() for f in patFragList])

        # Select proper fragments from mother/father
        newMatFragList, newPatFragList = [], []
        newMatFragAtms, newPatFragAtms = Atoms(), Atoms()
        for i in range(len(matDist)):
            if matDist[i] <= 0:
                matFrag = matFragList[i]
                newMatFragList.append(matFrag)
                newMatFragAtms.extend(matAds[matFrag])
        for i in range(len(patDist)):
            if patDist[i] > 0:
                patFrag = patFragList[i]
                newPatFragList.append(patFrag)
                newPatFragAtms.extend(patAds[patFrag])

        # Condense and sort each
        if len(newMatFragAtms) > 0:
            newMatFragList = frag.remake(newMatFragList,frag.flatsort(newMatFragList),range(len(newMatFragAtms))) 
            newMatFragAtms = sort(newMatFragAtms)
            #print([newMatFragAtms[f].get_chemical_formula() for f in newMatFragList])
        if len(newPatFragAtms) > 0:
            newPatFragList = frag.remake(newPatFragList,frag.flatsort(newPatFragList),range(len(newPatFragAtms))) 
            newPatFragAtms = sort(newPatFragAtms) 
            #print([newPatFragAtms[f].get_chemical_formula() for f in newPatFragList])

        # Fuse together
        newStructure = newMatFragList + newPatFragList
        newContents = frag.flatten(newMatFragList) + [len(frag.flatten(newMatFragList)) + i for f in newPatFragList for i in f]
        newFragList = frag.refill(newStructure,newContents)
        newFragAtms = newMatFragAtms + newPatFragAtms
        #print([newFragAtms[f].get_chemical_formula() for f in newFragList])

        # Transpose indices up to all interface atoms and set fragment atoms
        kid_fragList = frag.transposeUp(newFragList, childSurf)
        childSurf.set_fragList(kid_fragList)
        newFragAtms.info['adsorbate_fragments'] = frag.transposeUp(newFragList, childSurf)
        childSurf.set_adsAtoms_frag(newFragAtms)
        childSurf.wrap()

        # Evaluate quality of kid structure
        isBADSTRUCTURE = childSurf.has_badContact(tolerance=tolerance)

        # check bonds
        if not isBADSTRUCTURE:
            if bondRejList is not None:
                mySymb = childSurf.get_chemical_symbols()
                bps = geom.get_bondpairs(childSurf.get_allAtoms(), min(1-tolerance/2, 1-tolerance+0.2))
                #myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
                bps = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in bps]
                for br in bondRejList:
                    if br in bps or [br[1],br[0]] in bps:
                        #print('Rejected bond exists: ', br)
                        isBADSTRUCTURE = True
                        break
        #         if isBADSTRUCTURE:
        #             print('b', end='')
        # else:
        #     print(f'c', end='')

        if n_trial > int(200/tolerance):
            isBADSTRUCTURE = False
            childSurf = None
            print(f' FAIL ', end='')
            break
    if childSurf is not None:
        print('\nOffspring is created at attempt #%i\t|Tolerance = %.3f'%\
            (n_trial, tolerance))
    return childSurf
