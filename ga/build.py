

from gaia.utils import geom
from ase import Atoms
import numpy as np

def grow_adatom(
    interfc, addElemList,
    xLim=None, yLim=None,zLim=None,
    bldaSigma=0.1, bldaScale=1,
    ljopt=False, ljstepsize=0.01, ljnsteps=250
    ):
    numAds = len(addElemList)
    badStructure = True
    n_place = 0
    n_attempts = 0
    while badStructure:
        ind_curr = 0
        bufferList = interfc.get_bufferList()
        tmpInterfc = interfc.copy()
        while len(tmpInterfc) < len(interfc) + numAds:
            i = np.random.choice(bufferList)
            coord = [0,0,-1]
            while not geom.is_withinPosLim(coord, xLim, yLim, zLim):
                growVec = geom.rand_direction()
                growVec *= geom.BLDA(
                    addElemList[ind_curr],
                    tmpInterfc.get_chemical_symbols()[i],
                    sigma=bldaSigma, scale=bldaScale
                    )
                coord = tmpInterfc.get_positions()[i] + growVec
                n_place += 1
            tmpInterfc.merge_adsorbate(Atoms(addElemList[ind_curr], [coord]))
            bufferList.append(len(tmpInterfc)-1)
            ind_curr += 1
        if geom.chk_bondlength(tmpInterfc.get_allAtoms(), radTol=0.5):
            badStructure = False
        n_attempts += 1
    print(' %i\tplacements| %i\ttabula rasa'%(n_place, n_attempts - 1))
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc