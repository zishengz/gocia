import sys, os
import numpy as np
import ase.io as fio
from gocia.interface import Interface

def rand_cut_2d(atoms1, atoms2):
    '''
    Returns the center of the COMs of the parents
        and a random slope in x-y plane.
    '''
    com1 = atoms1.get_center_of_mass()
    com2 = atoms2.get_center_of_mass()
    direction = 2*np.random.rand(2)-1
    com_share = (com1 + com2)[:2]/2
    return com_share, direction

def crossover_cnsSurf_2d(surf1, surf2):
    dist2lin = lambda p1, p2, p3: np.cross(p2-p1, p1-p3)/np.linalg.norm(p2-p1)
    center, direction = rand_cut_2d(surf1.get_adsAtoms(), surf2.get_adsAtoms())
    matPos, patPos = surf1.get_pos(), surf2.get_pos()
    goodPos, badPos = [], []
    print(center, direction)
    matDist = [dist2lin(center, center+direction, i) for i in matPos[:, :2]]
    patDist = [dist2lin(center, center+direction, i) for i in patPos[:, :2]]
    print(matDist)
    print(len(surf1))
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
    print(sum([len(l) for l in goodPos]))
    print(sum([len(l) for l in badPos]))

    return np.array([r[0] for r in goodPos])

sub = fio.read(sys.argv[1])
s1  = fio.read(sys.argv[2])
s2  = fio.read(sys.argv[3])

surf1 = Interface(tags='surf 1', subAtoms=sub, allAtoms=s1)
surf2 = Interface(tags='surf 2', subAtoms=sub, allAtoms=s2)

surf1.print()
surf2.print()

newPos = crossover_cnsSurf_2d(surf1, surf2)
newSurf = surf1.copy()
newSurf.set_allPos(newPos)
newSurf.write('new.vasp')

