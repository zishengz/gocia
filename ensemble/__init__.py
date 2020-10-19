import numpy as np
from ase.db import connect

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
        print('   |-Isomer %i'%len(list(set(isomerLabel))), end='\t:')
        for i in range(nStates):
            if isomerLabel[i] != -1:
                continue
            if simMat[currInd][i] == 1:
                isomerLabel[i] = currInd
                print(i, end=',')
        currInd += 1

        print('')
    isoUniq = set(isomerLabel)
    eneUniq = [eneArray[i] for i in isoUniq]
    magUniq = [magArray[i] for i in isoUniq]
    print('-'*60)
    print('%i isomers found from %s samples.\tEfficiency = %s'%\
        (len(isoUniq), nStates, '{:.2%}'.format(len(isoUniq)/nStates)))
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