import sys, os
import numpy as np
import ase.io as fio
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
    center, direction = split_2d(surf1.get_adsAtoms(), surf2.get_adsAtoms())
    matPos, patPos = surf1.get_pos(), surf2.get_pos()
    goodPos, badPos = [], []
    matDist = [dist2lin(center, center+direction, i) for i in matPos[:, :2]]
    patDist = [dist2lin(center, center+direction, i) for i in patPos[:, :2]]
    isBADSTRUCTURE = True
    n_trial = 0
    while isBADSTRUCTURE:
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
    print(' Offspring is created at attempt #%i\t|Tolerance = %.3f'%\
        (n_trial, tolerance))
    return childSurf

sub = fio.read(sys.argv[1])
s1  = fio.read(sys.argv[2])
s2  = fio.read(sys.argv[3])

surf1 = Interface(tags='surf 1', subAtoms=sub, allAtoms=s1)
surf2 = Interface(tags='surf 2', subAtoms=sub, allAtoms=s2)

# surf1.print()
# surf2.print()
traj = []
for i in range(10):
    newSurf = crossover_snsSurf_2d(surf1, surf2, tolerance=0.75)
    newSurf.preopt_hooke(toler=0.1)
    traj.append(newSurf.get_allAtoms())
fio.write('gen.xyz', traj)

