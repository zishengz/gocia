import numpy as np
from gocia import geom

flatten = lambda l: [item for sublist in l for item in sublist]

def sortedAdsFGP(intfc, cutoff = 0.9):
    tmpAtoms = intfc.get_allAtoms()
    adsList = list(intfc.get_adsList()[:,1])
    bufList = list(intfc.get_bufferList())
    adsBl = []
    for ads in adsList:
        neib = geom.get_neighbors(
            tmpAtoms,
            ads,
            bufList+adsList,
            scale = cutoff
            )
        if len(neib[0]) > 0:
            adsBl += list(tmpAtoms.get_distances(ads, neib[0], mic=True))
        adsBl.sort()
    return np.array(adsBl)

def sortedBufFGP(intfc, cutoff = 0.9):
    tmpAtoms = intfc.get_allAtoms()
#    adsList = list(intfc.get_adsList()[:,1])
    bufList = list(intfc.get_bufferList())
    bufBl = []
    for buf in bufList:
        neib = geom.get_neighbors(
            tmpAtoms,
            buf,
            bufList,#+adsList,
            scale = cutoff
            )
        if len(neib[0]) > 0:
            bufBl += list(tmpAtoms.get_distances(buf, neib[0], mic=True))
        bufBl.sort()
    return np.array(bufBl)

def sortedAdsBufFGP(intfc, cutoff = 0.9):
    return get_sortedAdsList(intfc, cutoff),\
           get_sortedBufList(intfc, cutoff)

def coordinationFGPv1(atoms, scale=1.25):
    '''
    Returns a 1x10 D array
    '''
    cn, nb, bl = geom.get_coordStatus(atoms, scale=scale)
    bl = np.array(flatten(bl))
    pos = atoms.get_positions()
    # X should only be enabled in grand canonical ensemble
    # or it would cause normalization problem
    myFGP = np.array([
        # atoms.get_masses().sum(),       # X:  total mass
        # len(atoms),                     # X:  number of atoms
        pos[:,2].mean(),                # 1:  mean of z positions
        pos[:,2].std(),                 # 2:  stdev of z positions
        pos[:,2].max()-pos[:,2].min(),  # 3:  slab thickness
        atoms.get_center_of_mass()[2],  # 4:  center of mass
        cn.mean(),                      # 5:  mean of CN
        cn.std(),                       # 6:  stdev of CN
        bl.mean(),                      # 7:  mean of bond length
        bl.std(),                       # 8:  stdev of bond length
        bl.min(),                       # 9: max of bond length
        bl.max(),                       # 10: max of bond length
    ])
    return myFGP
