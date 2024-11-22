import gocia.utils.linalg as la
import gocia.geom.fingerprint as fgp
import numpy as np

def get_sim_srtDist(a1, a2):
    # adapted from 10.1063/1.4886337
    if len(a1) != len(a2):
        print('they two sytems must be of the same size (#atoms)!')
        exit()
    p1 = a1.get_all_distances(mic=True).flatten()
    p2 = a2.get_all_distances(mic=True).flatten()
    p1, p2 = np.sort(p1), np.sort(p2)
    cum_diff = np.abs(p1 - p2)
    total_cum_diff = cum_diff.sum() / ((p1 + p2).sum()/2)
    max_diff = cum_diff.max()
    return total_cum_diff, max_diff

def srtDist_similar_zz(a1, a2, delta_rel=5e-4, d_max=0.25):
    if len(a1) != len(a2):
        return False
    else:
        total_cum_diff, max_diff = get_sim_srtDist(a1, a2)
        
        if total_cum_diff < delta_rel and max_diff < d_max:
            return True
        else:
            return False


def comp_srtDist_zz(a1, traj, delta_rel=8e-4, d_max=0.4):
    simList = [0]*len(traj)
    for i in range(len(traj)):
        if srtDist_similar_zz(a1, traj[i], delta_rel=delta_rel, d_max=d_max):
            simList[i] = 1
    return simList


def compAll_srtDist_zz(traj, delta_rel=8e-4, d_max=0.4):
    simMat = []
    for i in range(len(traj)):
        print('checking isomer %i' % i)
        simMat.append(comp_srtDist_zz(
            traj[i], traj, delta_rel=delta_rel, d_max=d_max))
    return np.array(simMat)

# Above can do one-by-one and one-by-many checks
# Below functions only apply for a whole ensemble


def num_passed(simMat):
    return (len(simMat[simMat == 1]) - len(simMat))/2


def compare_geom(traj, cutoff=0.1,
                 myFGP=fgp.coordinationFGPv1):
    '''
    cutoff > 0.9 is usually good
    '''
    print(' * Analyzing similarity in GEOM. FINGERPRINT...\tCutoff = %.3f' % (cutoff))
    fgpArr = np.array([myFGP(a) for a in traj])
    geomDiff = la.cosSimMatrix(fgpArr)
    geomDiff[geomDiff <= cutoff] = -1.0
    geomDiff[geomDiff > cutoff] = 0
    geomDiff[geomDiff == -1.0] = 1
    print('   |- Pairs passed: %i' % num_passed(geomDiff))
    return geomDiff


def compare_ene(ene, cutoff=0.1):
    '''
    return value of 1 means pass
    '''
    print(' * Analyzing similarity in ELECTRONIC ENERGY...\tCutoff = %.3f' % (cutoff))
    if type(ene) is list:
        ene = np.array(ene)
    eneDiff = np.abs(la.diffMatrix(ene))
    eneDiff[eneDiff <= cutoff] = -1.0
    eneDiff[eneDiff > cutoff] = 0
    eneDiff[eneDiff == -1.0] = 1
    print('   |- Pairs passed: %i' % (num_passed(eneDiff)))
    return eneDiff


def compare_posEig(traj, cutoff):
    print(' * Analyzing similarity in DISTANCE MATRIX  ...\tCutoff = %.3f' % (cutoff))
    allEig = np.array([fgp.posMatEigenFGP(i) for i in traj])
    allEigDist = np.array([la.euclDist(allEig, e, axis=1) for e in allEig])
    allEigDist[allEigDist <= 1] = -1
    allEigDist[allEigDist > 1] = 0
    allEigDist[allEigDist == -1] = 1
    print('   |- Pairs passed: %i' % (num_passed(allEigDist)))
    return allEigDist


def bothSim(simMat1, simMat2):
    return la.normalize_mat(simMat1 + simMat2)


def compare_dual(traj, geomCutoff=0.1, enerCutoff=0.1,
                 myFGP=fgp.coordinationFGPv1):
    ene = [a.info['eV'] for a in traj]
    simGeom = compare_geom(traj, geomCutoff)
    simEner = compare_ene(ene, enerCutoff)
    return bothSim(simEner, simGeom)
