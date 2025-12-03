import numpy as np
from ase.db import connect


def clusterIsomer(
    trajAtoms,
    simMat,
    outName='sort'
):
    print(' * Detecting unique isomers...')
    nStates = len(simMat)
    isomerLabel = [-1]*nStates
    currIso, currInd = [], 0
    while -1 in isomerLabel:
        if currInd >= nStates:
            break
        if isomerLabel[currInd] != -1:
            currInd += 1
            continue
        print('   |-Isomer %i' % len(list(set(isomerLabel))), end='\t:')
        for i in range(nStates):
            if isomerLabel[i] != -1:
                continue
            if simMat[currInd][i] == 1:
                isomerLabel[i] = currInd
                print(i, end=',')
        currInd += 1
        print('')

    # Handle the bad or other-format outputs
    if 'mag' in trajAtoms[0].info.keys():
        magArray = [a.info['mag'] for a in trajAtoms]
    else:
        magArray = [0]*len(trajAtoms)
    if 'eV' in trajAtoms[0].info.keys():
        eneArray = [a.info['eV'] for a in trajAtoms]
    else:
        eneArray = [a.get_potential_energy() for a in trajAtoms]

    isoUniq = set(isomerLabel)
    eneUniq = [eneArray[i] for i in isoUniq]
    magUniq = [magArray[i] for i in isoUniq]
    print('-'*76)
    print('   %i\tisomers found from %s\tsamples. Oversampling = %s' %
          (len(isoUniq), nStates, '{:.2%}'.format(nStates/len(isoUniq)-1)))
    sortedTraj = []
    for v in sorted(eneUniq):
        sortedTraj.append(trajAtoms[list(eneArray).index(v)])

    with connect('%s.db' % outName, append=False) as srtDB:
        for i in range(len(sortedTraj)):
            srtDB.write(
                sortedTraj[i],
                done=1,
                eV=sorted(eneUniq)[i],
                mag=magUniq[eneUniq.index(sorted(eneUniq)[i])],
                eV2GM=sorted(eneUniq)[i] - min(eneUniq)
            )


def clusterIsomerAbs(
    trajAtoms,  # Must containg energy information
    eneToler=0.05,
    geomToler1=5e-4,
    geomToler2=0.5,
    outName='sort'
):
    import gocia.ensemble.comparator as comp

    # Handle the bad or other-format outputs
    if 'mag' in trajAtoms[0].info.keys():
        magArray = [a.info['mag'] for a in trajAtoms]
    else:
        magArray = [0]*len(trajAtoms)
    if 'eV' in trajAtoms[0].info.keys():
        eneArray = [a.info['eV'] for a in trajAtoms]
    else:
        eneArray = [a.get_potential_energy() for a in trajAtoms]
    if 'adsFrags' in trajAtoms[0].info.keys():
        fragArray = [a.info['adsFrags'] for a in trajAtoms]
    else:
        fragArray = []
        print('WARNING: No frag list is found in the databse')
        #print('Uh oh, gotta get frag list somehow into gocia.ensemble.clusterIsomerAbs')

    print(' * Detecting unique isomers...')
    nStates = len(trajAtoms)
    isomerLabel = [-1]*nStates
    currIso, currInd = [], 0
    while -1 in isomerLabel:
        if currInd >= nStates:
            break
        if isomerLabel[currInd] != -1:
            currInd += 1
            continue
        print('   |-Isomer %i' % len(list(set(isomerLabel))), end='\t:')
        for i in range(nStates):
            if isomerLabel[i] != -1:
                continue
            # if simMat[currInd][i] == 1:
            #     isomerLabel[i] = currInd
            if abs(eneArray[i] - eneArray[currInd]) < eneToler and\
                comp.srtDist_similar_zz(trajAtoms[i], trajAtoms[currInd], delta_rel=geomToler1, d_max=geomToler2) and\
                    trajAtoms[i].get_chemical_formula() == trajAtoms[currInd].get_chemical_formula():
                isomerLabel[i] = currInd
                print(f'{i} ({trajAtoms[i].get_chemical_formula()})', end=', ')
        currInd += 1
        print('')

    isoUniq = set(isomerLabel)
    eneUniq = [eneArray[i] for i in isoUniq]
    magUniq = [magArray[i] for i in isoUniq]
    try:
        fragUniq = [fragArray[i] for i in isoUniq]
    except:
        pass
    countUniq = [isomerLabel.count(i) for i in isoUniq]
#    print(countUniq)
    print('-'*76)
    print('   %i\tisomers found from %s\tsamples. Oversampling = %s' %
          (len(isoUniq), nStates, '{:.2%}'.format(nStates/len(isoUniq)-1)))
    sortedTraj = []
    for v in sorted(eneUniq):
        sortedTraj.append(trajAtoms[list(eneArray).index(v)])

    with connect('%s.db' % outName, append=False) as srtDB:
        for i in range(len(sortedTraj)):
            srtDB.write(
                sortedTraj[i],
                done=1,
                eV=sorted(eneUniq)[i],
                eV2GM=sorted(eneUniq)[i] - min(eneUniq),
                mag=magUniq[eneUniq.index(sorted(eneUniq)[i])],
                counts=countUniq[eneUniq.index(sorted(eneUniq)[i])],
#                adsFrags=fragUniq[eneUniq.index(sorted(eneUniq)[i])],
            )


def clusterIsomerAbs_opt(
    trajAtoms,  # Must containg energy information
    eneToler=0.05,
    geomToler1=5e-4,
    geomToler2=0.5,
    outName='sort'
):
    import gocia.ensemble.comparator as comp

    # Handle the bad or other-format outputs
    if 'mag' in trajAtoms[0].info.keys():
        magArray = [a.info['mag'] for a in trajAtoms]
    else:
        magArray = [0]*len(trajAtoms)
    if 'eV' in trajAtoms[0].info.keys():
        eneArray = [a.info['eV'] for a in trajAtoms]
    else:
        eneArray = [a.get_potential_energy() for a in trajAtoms]
    if 'adsFrags' in trajAtoms[0].info.keys():
        fragArray = [a.info['adsFrags'] for a in trajAtoms]
    else:
        fragArray = []
        print('WARNING: No frag list is found in the databse')
        #print('Uh oh, gotta get frag list somehow into gocia.ensemble.clusterIsomerAbs')

    print(' * Computing all similarity descriptors...')
    list_sortdist = []
    for s in trajAtoms:
        tmp = s.get_all_distances(mic=True).flatten()
        list_sortdist.append(np.sort(tmp))
    del tmp

    print(' * Detecting unique isomers...')
    nStates = len(trajAtoms)
    isomerLabel = [-1]*nStates
    currIso, currInd = [], 0
    while -1 in isomerLabel:
        if currInd >= nStates:
            break
        if isomerLabel[currInd] != -1:
            currInd += 1
            continue
        print('   |-Isomer %i' % len(list(set(isomerLabel))), end='\t:')
        for i in range(nStates):
            if isomerLabel[i] != -1:
                continue
            if trajAtoms[i].get_chemical_formula() != trajAtoms[currInd].get_chemical_formula():
                continue
            if abs(eneArray[i] - eneArray[currInd]) > eneToler:
                continue
            p1 = list_sortdist[i]
            p2 = list_sortdist[currInd]
            p1, p2 = np.sort(p1), np.sort(p2)
            cum_diff = np.abs(p1 - p2)
            total_cum_diff = cum_diff.sum() / ((p1.sum() + p2.sum())/2)
            max_diff = cum_diff.max()
            if total_cum_diff < geomToler1 and max_diff < geomToler2:
                isomerLabel[i] = currInd
                print(f'{i} ({trajAtoms[i].get_chemical_formula()})', end=', ')
        currInd += 1
        print('')

    isoUniq = set(isomerLabel)
    eneUniq = [eneArray[i] for i in isoUniq]
    magUniq = [magArray[i] for i in isoUniq]
    try:
        fragUniq = [fragArray[i] for i in isoUniq]
    except:
        pass
    countUniq = [isomerLabel.count(i) for i in isoUniq]
#    print(countUniq)
    print('-'*76)
    print('   %i\tisomers found from %s\tsamples. Oversampling = %s' %
          (len(isoUniq), nStates, '{:.2%}'.format(nStates/len(isoUniq)-1)))
    sortedTraj = []
    for v in sorted(eneUniq):
        sortedTraj.append(trajAtoms[list(eneArray).index(v)])

    with connect('%s.db' % outName, append=False) as srtDB:
        for i in range(len(sortedTraj)):
            srtDB.write(
                sortedTraj[i],
                done=1,
                eV=sorted(eneUniq)[i],
                eV2GM=sorted(eneUniq)[i] - min(eneUniq),
                mag=magUniq[eneUniq.index(sorted(eneUniq)[i])],
                counts=countUniq[eneUniq.index(sorted(eneUniq)[i])],
#                adsFrags=fragUniq[eneUniq.index(sorted(eneUniq)[i])],
            )


def clusterIsomerAbs_GC(
    trajAtoms,  # Must containg energy information
    eneToler=0.05,
    geomToler1=5e-4,
    geomToler2=0.25,
    outName='sort'
):
    import gocia.ensemble.comparator as comp

    # Handle the bad or other-format outputs
    if 'mag' in trajAtoms[0].info.keys():
        magArray = [a.info['mag'] for a in trajAtoms]
    else:
        magArray = [0]*len(trajAtoms)
    if 'grandPot' in trajAtoms[0].info.keys():
        eneArray = [a.info['grandPot'] for a in trajAtoms]
    else:
        try:
            eneArray = [a.get_potential_energy() for a in trajAtoms]
        except:
            eneArray = [a.info['eV'] for a in trajAtoms]
    rawEneArray = [a.info['eV'] for a in trajAtoms]

    print(' * Detecting unique isomers...')
    nStates = len(trajAtoms)
    isomerLabel = [-1]*nStates
    currIso, currInd = [], 0
    while -1 in isomerLabel:
        if currInd >= nStates:
            break
        if isomerLabel[currInd] != -1:
            currInd += 1
            continue
        print('   |-Isomer %i' % len(list(set(isomerLabel))), end='\t:')
        for i in range(nStates):
            if isomerLabel[i] != -1:
                continue
            # if simMat[currInd][i] == 1:
            #     isomerLabel[i] = currInd
            if abs(eneArray[i] - eneArray[currInd]) < eneToler and\
                    comp.srtDist_similar_zz(trajAtoms[i], trajAtoms[currInd],
                                            delta_rel=geomToler1, d_max=geomToler2) and\
                    trajAtoms[i].get_chemical_formula() == trajAtoms[currInd].get_chemical_formula():
                isomerLabel[i] = currInd
                print(f'{i} ({trajAtoms[i].get_chemical_formula()})', end=', ')
        currInd += 1
        print('')

    isoUniq = set(isomerLabel)
    eneUniq = [eneArray[i] for i in isoUniq]
    magUniq = [magArray[i] for i in isoUniq]
    rawUniq = [rawEneArray[i] for i in isoUniq]
    countUniq = [isomerLabel.count(i) for i in isoUniq]
#    print(countUniq)
    print('-'*76)
    print('   %i\tisomers found from %s\tsamples. Oversampling = %s' %
          (len(isoUniq), nStates, '{:.2%}'.format(nStates/len(isoUniq)-1)))
    sortedTraj = []
    for v in sorted(eneUniq):
        sortedTraj.append(trajAtoms[list(eneArray).index(v)])

    with connect('%s.db' % outName, append=False) as srtDB:
        for i in range(len(sortedTraj)):
            srtDB.write(
                sortedTraj[i],
                done=1,
                grandPot=sorted(eneUniq)[i],
                eV2GM=sorted(eneUniq)[i] - min(eneUniq),
                eV=rawUniq[eneUniq.index(sorted(eneUniq)[i])],
                mag=magUniq[eneUniq.index(sorted(eneUniq)[i])],
                counts=countUniq[eneUniq.index(sorted(eneUniq)[i])],
            )
