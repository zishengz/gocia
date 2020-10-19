
import sys, os
import numpy as np
import ase.io as fio
from ase.db import connect
import gocia.geom.fingerprint as fgp

def get_property(dbobj, keyVal):
    '''
    Returns 1-D list
    '''
    tmp = []
    for r in dbobj.select():
        tmp.append(r[keyVal])
    return tmp

def get_fgp(dbobj, myfgp):
    tmp = []
    for r in dbobj.select():
        tmp.append(myfgp(r.toatoms()))
    return np.array(tmp)

def normalize_mat(matrix):
    return (matrix - matrix.min())\
        / (matrix.max() - matrix.min())

def normalize_fgp(fgpArray):
    return (fgpArray - fgpArray.min(axis=0))\
        / (fgpArray.max(axis=0) - fgpArray.min(axis=0))

def cosSimMatrix(fgpArray):
    from sklearn.metrics.pairwise import cosine_similarity
    if fgpArray.max() == 1 and fgpArray.min() == 0:
        return cosine_similarity(fgpArray)
    else:
        return cosine_similarity(normalize_fgp(fgpArray))

def diffMatrix(d1Array):
    return d1Array[:, np.newaxis] - d1Array

def heatmap(matrix, baseName):
    import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
    dimension = len(matrix)
    fig, ax = plt.subplots()
    im = ax.imshow(matrix, cmap=cm.Blues)
    fig.colorbar(im)
    ax.tick_params(top=True, bottom=False,
                    labeltop=True, labelbottom=False)
    ax.set_xticks(np.arange(dimension))
    ax.set_yticks(np.arange(dimension))
    # ... and label them with the respective list entries
    ax.set_xticklabels([str(i) for i in np.arange(dimension)], fontsize='8')
    ax.set_yticklabels([str(i) for i in np.arange(dimension)], fontsize='8')
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right",
            rotation_mode="anchor")
    fig.tight_layout()
    plt.savefig(
        '%s_heat.png'%(baseName),
        dpi=200
        )

def clusterIsomer(
    trajAtoms,
    simMat,
    eneArray,
    magArray,
    outName='sort'
    ):
    nStates = len(simMat)
    isomerLabel = [-1]*nStates
    currIso, currInd = [], 0
    while -1 in isomerLabel:
        if currInd >= nStates:
            break
        if isomerLabel[currInd] != -1:
            currInd += 1
            continue
        print('Isomers of %i'%currInd, end=':\t')
        for i in range(nStates):
            if isomerLabel[i] != -1:
                continue
            if simMat[currInd][i] == 1:
                isomerLabel[i] = currInd
                print(i, end=',')
        currInd += 1
        print('')
    print(isomerLabel)
    isoUniq = set(isomerLabel)
    eneUniq = [eneArray[i] for i in isoUniq]
    magUniq = [magArray[i] for i in isoUniq]
    print(eneUniq)
    print('Number of isomers: %i'%len(isoUniq))
    sortedTraj = []
    for v in sorted(eneUniq):
        sortedTraj.append(trajAtoms[list(eneArray).index(v)])

    srtDB = connect('%s.db'%outName, append=False)
    for i in range(len(sortedTraj)):
        srtDB.write(
            sortedTraj[i],
            eV = sorted(eneUniq)[i],
            mag = magUniq[eneUniq.index(sorted(eneUniq)[i])],
            eV2GM = sorted(eneUniq)[i] - min(eneUniq)
        )

# $ python test.py raw-N1.db
# # Timing: 5.7539 s
trajDBName = sys.argv[1]
geomCutoff = 0.9
enerCutoff = 0.1

print(' |- Extracting data from database: %s'%trajDBName)
rawDB = connect(trajDBName)
fgpArr = get_fgp(rawDB, fgp.coordinationFGPv1)
ene = np.array(get_property(rawDB, 'eV'))
mag = np.array(get_property(rawDB, 'mag'))
traj = fio.read(trajDBName, index=':')

print(' |- Analyzing geometry similarity')
geomSim = cosSimMatrix(fgpArr)
geomSim[geomSim <  geomCutoff] = 0
geomSim[geomSim >= geomCutoff] = 1.0
heatmap(geomSim, 'geom'+str(geomCutoff))

print(' |- Analyzing energy similarity')
eneDiff = np.abs(diffMatrix(ene))
eneDiff[eneDiff <= enerCutoff] = -1.0
eneDiff[eneDiff >  enerCutoff] = 0
eneDiff[eneDiff == -1.0] = 1
heatmap(eneDiff, 'ener'+str(enerCutoff))

print(' |- Analyzing overall similarity')
bothPass = normalize_mat(eneDiff + geomSim)
heatmap(bothPass, 'both_G%sE%s'%\
    (str(geomCutoff), str(enerCutoff)))

print(' |- Cluster ensemble into isomers')
clusterIsomer(
    traj,
    bothPass,
    ene,
    mag,
)


