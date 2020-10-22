

from gocia import geom
from ase import Atoms
import numpy as np
import random

def grow_adatom(
    interfc, addElemList,
    xLim=None, yLim=None,zLim=None,
    sampZEnhance=None,
    bldaSigma=0.1, bldaScale=1, toler=0.5,
    doShuffle=False,
    rattle=False, rattleStdev = 0.05,rattleZEnhance=False,
    sameElemPenalty = 0,
    cnCount=False,cnToler = 0.5,
    ljopt=False, ljstepsize=0.01, ljnsteps=400
    ):
    numAds = len(addElemList)
    badStructure = True
    n_place = 0
    n_attempts = 0
    while badStructure:
        if doShuffle:
            random.shuffle(addElemList)
        ind_curr = 0
        tmpInterfc = interfc.copy()
        if rattle:
            tmpInterfc.rattle(rattleStdev, zEnhance=rattleZEnhance)
        while len(tmpInterfc) < len(interfc) + numAds and ind_curr < len(addElemList):
            optList = tmpInterfc.get_optList()
            weights = np.ones(len(optList))
            # Same-element penalty
            # value of 0 has no effect on weighting
            # set to >0 to avoid same elem sampling:
            #       1 -> same 33%, diff 66%
            #       3 -> same 25%, diff 75%
            # set to -1 for pure cluster growth (naive implementation)
            optElem = [tmpInterfc.get_chemical_symbols()[i] for i in optList]
            mult = np.ones(len(optList)) +\
                   np.array([sameElemPenalty if s != addElemList[ind_curr] else 0 for s in optElem])
            if np.count_nonzero(mult) == 0:
                mult = np.ones(len(optList))
            weights *= mult
            # Z-weighted sampling enhancement
            if sampZEnhance is not None:
                # sampAEnhance = 1 makes p(zmax) = p(zmin)*2
                optZ = tmpInterfc.get_pos()[:,2][[i for i in optList]]
                if optZ.max() - optZ.min() != 0:
                    mult = 1 + (optZ - optZ.min())/(optZ.max() - optZ.min()) * sampZEnhance
                    weights *= mult
            weights /= weights.sum()
            i = np.random.choice(optList, p=weights)
            # Coordination-adapted sampling
            if cnCount:
                cn = geom.get_coordStatus(tmpInterfc.get_allAtoms())[0]
                if cn.max() - cn.min() != 0:
                    cn = (cn[i] - cn.min()) / (cn.max() - cn.min())
                    if np.random.rand() < cn * cnToler: continue
#             if addElemList[ind_curr] == tmpInterfc.get_chemical_symbols()[i]:
#                 if np.random.rand() >  1/2 + sameElemPenalty: continue
#             else:
#                 if np.random.rand() <= 1/2 + sameElemPenalty: continue
# #                if np.random.rand() < sameElemPenalty: continue
            coord = [0,0,-1000]
            while not geom.is_withinPosLim(coord, xLim, yLim, zLim):
                growVec = geom.rand_direction()
                blda = geom.BLDA(
                    addElemList[ind_curr],
                    tmpInterfc.get_chemical_symbols()[i],
                    sigma=bldaSigma, scale=bldaScale
                    )
                growVec *= blda
                coord = tmpInterfc.get_pos()[i] + growVec
#                print(blda, growVec, coord)
#                print(i, tmpInterfc.get_pos()[i])
                n_place += 1
#            print('pass')
            tmpInterfc.merge_adsorbate(Atoms(addElemList[ind_curr], [coord]))
            ind_curr += 1
        if tmpInterfc.has_badContact(tolerance = toler):
            badStructure = False
        n_attempts += 1
    tmpInterfc.sort()
    print('%i\tplacements| %i\ttabula rasa'%(n_place, n_attempts - 1))
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc
