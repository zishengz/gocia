import gocia.utils.linalg as la
import gocia.geom.fingerprint as fgp
import numpy as np

def num_passed(simMat):
    return (len(simMat[simMat==1]) - len(simMat))/2

def compare_geom(traj, cutoff=0.9,
                 myFGP=fgp.coordinationFGPv1):
    '''
    cutoff > 0.9 is usually good
    '''
    print(' * Analyzing similarity in GEOMETRY...\tCutoff = %.3f'%(cutoff))
    fgpArr = np.array([myFGP(a) for a in traj])
    geomSim = la.cosSimMatrix(fgpArr)
    geomSim[geomSim <  cutoff] = 0
    geomSim[geomSim >= cutoff] = 1.0
    print('   |- Pairs passed: %i'%num_passed(geomSim))
    return geomSim

def compare_ene(ene,cutoff=0.02):
    '''
    return value of 1 means pass
    '''
    print(' * Analyzing similarity in ENERGY...\tCutoff = %.3f'%(cutoff))
    if type(ene) is list:
        ene = np.array(ene)
    eneDiff = np.abs(la.diffMatrix(ene))
    eneDiff[eneDiff <= cutoff] = -1.0
    eneDiff[eneDiff >  cutoff] = 0
    eneDiff[eneDiff == -1.0] = 1
    print('   |- Pairs passed: %i'%(num_passed(eneDiff)))
    return eneDiff

def bothSim(simMat1, simMat2):
    return la.normalize_mat(simMat1 + simMat2)

def compare_dual(traj, geomCutoff=0.9, enerCutoff=0.02,
                 myFGP=fgp.coordinationFGPv1):
    ene = [a.info['eV'] for a in traj]
    simGeom = compare_geom(traj, geomCutoff)
    simEner = compare_ene(ene, enerCutoff)
    return bothSim(simEner, simGeom)