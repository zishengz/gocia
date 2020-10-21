import gocia.utils.linalg as la
import gocia.geom.fingerprint as fgp
import numpy as np

def num_passed(simMat):
    return (len(simMat[simMat==1]) - len(simMat))/2

def compare_geom(traj, cutoff=0.1, 
                 myFGP=fgp.coordinationFGPv1):
    '''
    cutoff > 0.9 is usually good
    '''
    print(' * Analyzing similarity in GEOM. FINGERPRINT...\tCutoff = %.3f'%(cutoff))
    fgpArr = np.array([myFGP(a) for a in traj])
    geomDiff = la.cosSimMatrix(fgpArr)
    geomDiff[geomDiff <=  cutoff] = -1.0
    geomDiff[geomDiff > cutoff] = 0
    geomDiff[geomDiff == -1.0] = 1
    print('   |- Pairs passed: %i'%num_passed(geomDiff))
    return geomDiff

def compare_ene(ene,cutoff=0.1):
    '''
    return value of 1 means pass
    '''
    print(' * Analyzing similarity in ELECTRONIC ENERGY...\tCutoff = %.3f'%(cutoff))
    if type(ene) is list:
        ene = np.array(ene)
    eneDiff = np.abs(la.diffMatrix(ene))
    eneDiff[eneDiff <= cutoff] = -1.0
    eneDiff[eneDiff >  cutoff] = 0
    eneDiff[eneDiff == -1.0] = 1
    print('   |- Pairs passed: %i'%(num_passed(eneDiff)))
    return eneDiff

def compare_posEig(traj, cutoff):
    print(' * Analyzing similarity in DISTANCE MATRIX  ...\tCutoff = %.3f'%(cutoff))
    allEig = np.array([fgp.posMatEigenFGP(i) for i in traj])
    allEigDist = np.array([la.euclDist(allEig, e, axis=1) for e in allEig])
    allEigDist[allEigDist<= 1] = -1
    allEigDist[allEigDist> 1] = 0
    allEigDist[allEigDist == -1] = 1
    print('   |- Pairs passed: %i'%(num_passed(allEigDist)))
    return allEigDist

def bothSim(simMat1, simMat2):
    return la.normalize_mat(simMat1 + simMat2)

def compare_dual(traj, geomCutoff=0.1, enerCutoff=0.1,
                 myFGP=fgp.coordinationFGPv1):
    ene = [a.info['eV'] for a in traj]
    simGeom = compare_geom(traj, geomCutoff)
    simEner = compare_ene(ene, enerCutoff)
    return bothSim(simEner, simGeom)