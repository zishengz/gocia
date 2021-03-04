
from ase.io import read, write
from ase.data import covalent_radii
from scipy.spatial.transform import Rotation as R
import numpy as np
import random
import sys
from gocia import geom

def growAdsorbate(clusterAtoms, adsorbateList, region_radii):
    com = clusterAtoms.get_center_of_mass()
    tmpAtoms = clusterAtoms.copy()
    # Shuffle the order of added adsorbates
    random.shuffle(adsorbateList)
    ads_to_add = adsorbateList.copy()
    while len(ads_to_add) > 0:
        tmpAds = ads_to_add[0].copy()
        del tmpAds[-1]
        ads_site_index = np.random.choice(range(len(clusterAtoms)))
        ads_site_pos = clusterAtoms.get_positions()[ads_site_index]
        vec_grow = ads_site_pos - com
        dist = np.linalg.norm(vec_grow)
        vec_grow /= dist
        rot_axis = np.cross(vec_grow, np.array([0,0,1]))
        rot_axis /= np.linalg.norm(rot_axis)
        rot_angle = np.arccos(np.dot(vec_grow, np.array([0,0,1])))
        rot = R.from_rotvec(- rot_angle * rot_axis)
        pos_ads = tmpAds.get_positions()
        pos_ads = rot.apply(pos_ads)
        pos_ads += dist*vec_grow
        # TODO
        # random rotation of adsorbate
        # adjusting distance
        tmpAds.set_positions(pos_ads)
        tmpAtoms += tmpAds
        break
    return tmpAtoms

def get_covalRadii(atoms):
    return [covalent_radii[i] for i in atoms.numbers]

def get_contactMat(atoms, scale=1.0):
    '''
    returns the standard covalent bondlength between all atoms
    '''
    cvRad = get_covalRadii(atoms)
    contact = np.tile(cvRad, [len(cvRad), 1])
    contact += contact.T
    contact *= scale
    np.fill_diagonal(contact, 0)
    return contact

def count_badContact(atoms, tolerance=0):
    diff = atoms.get_all_distances(mic=True) - get_contactMat(atoms, scale=1-tolerance)
    #diff = atoms.get_all_distances(mic=True) - get_contactMat(atoms, scale=1) + tolerance
    diff[diff>0] = 0
    diff[diff<0] = 1
    return diff.sum()

def growAdsorbate_landing(clusterAtoms, adsorbateList, region_radii=[1,1,1], toler=0.1):
    com = clusterAtoms.get_center_of_mass()
    dist2com = np.sqrt(((clusterAtoms.get_positions()-com)**2).sum(axis=1))
    tmpAtoms = clusterAtoms.copy()
    # Shuffle the order of added adsorbates
    random.shuffle(adsorbateList)
    ads_to_add = adsorbateList.copy()
    n_added = 0
    n_trials = 0
    while n_added < len(ads_to_add):
        n_trials += 1
#        print('Adsorbates to add: %i/%i'%(n_added, len(ads_to_add)))
        tmpAds = ads_to_add[0].copy()
        del tmpAds[-1]
        vec_grow = geom.rand_direction() * np.array(region_radii)
        vec_grow /= np.linalg.norm(vec_grow)
        rot_axis = np.cross(vec_grow, np.array([0,0,1]))
        rot_axis /= np.linalg.norm(rot_axis)
        rot_angle = np.arccos(np.dot(vec_grow, np.array([0,0,1])))
        rot = R.from_rotvec(- rot_angle * rot_axis)
        pos_ads = tmpAds.get_positions()
        pos_ads = rot.apply(pos_ads)
        dist = max(dist2com) + 3
        badStructure = True
        trj=[]
        myContacts=count_badContact(tmpAtoms, toler)
        while badStructure:
            pos_ads_test = pos_ads + dist*vec_grow
            tmpClus = clusterAtoms.copy()
            tmpAtoms_test = tmpAtoms.copy() 
            tmpAds_test = tmpAds.copy()
            tmpAds_test.set_positions(pos_ads_test)
            tmpAtoms_test += tmpAds_test
            tmpClus += tmpAds_test
            trj.append(tmpAtoms_test)
            clusContacts = count_badContact(tmpClus, toler)
            if clusContacts > 0:
                if n_added == 0:
                    badStructure = False
                    break
#                print(myContacts, clusContacts, count_badContact(tmpAtoms_test, toler))
                if count_badContact(tmpAtoms_test, toler) ==clusContacts+myContacts:
                    badStructure = False
                    break
                else:
                    break
            dist -= 0.025
            if dist < min(dist2com):
                break
        if not badStructure:
            tmpAtoms = tmpAtoms_test.copy()
            n_added += 1
        write('test.xyz', trj)
    print('Acceptance rate: %.2f%%'%(len(adsorbateList)/n_trials*100))
    return tmpAtoms

clus = read('./Pt13-gm.xyz')
ads = read('./H.xyz')

traj = []
for i in range(10):
    print(i, end='\t')
    clus_ads = growAdsorbate_landing(clus, [ads]*8, [1,1,1], 0.1)
    traj.append(clus_ads)
write('ini.xyz', traj)



