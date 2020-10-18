import numpy as np
from gocia import geom

def get_sortedAdsList(intfc, cutoff = 0.9):
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

def get_sortedBufList(intfc, cutoff = 0.9):
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

def get_sortedBondLists(intfc, cutoff = 0.9):
    return get_sortedAdsList(intfc, cutoff),\
           get_sortedBufList(intfc, cutoff)
