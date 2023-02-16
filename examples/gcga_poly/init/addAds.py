from gocia.interface import Interface
from gocia import geom
from ase.data import vdw_radii, covalent_radii
from ase.io import read, write
import numpy as np
import random
from scipy.spatial.transform import Rotation as R
from ase.db import connect

###
#   addAds.py
# 
#   This script adds adsorbates to a substrate, tracks the atoms that compose each added adsorbate, 
#   and generates a database with the new structures and lists of atom indices in each adsorbate fragment.
#
#   Modified by Winston Gee from solvSamp.py created by Zisheng Zhang
###

# Helper functions
def get_vdwRadii(atoms):
    return [vdw_radii[i] for i in atoms.numbers]

def get_covRadii(atoms):
    return [covalent_radii[i] for i in atoms.numbers]

def get_contactMat(atoms, scale=1.0):
    '''
    returns the standard covalent bondlength between all atoms
    '''
    myRad = get_covRadii(atoms)
    contact = np.tile(myRad, [len(myRad), 1])
    contact += contact.T
    contact *= scale
    np.fill_diagonal(contact, 0)
    return contact

def count_badContact(atoms, tolerance=-0.5):
    #diff = atoms.get_all_distances(mic=True) - get_contactMat(atoms, scale=1-tolerance)
    diff = atoms.get_all_distances(mic=True) - get_contactMat(atoms, scale=1)
    count_tooClose = np.count_nonzero(diff < tolerance)/2
    return count_tooClose

def count_bondContact(atoms, tolerLims=[-0.3, 0.2]):
    #diff = atoms.get_all_distances(mic=True) - get_contactMat(atoms, scale=1-tolerance)
    diff = atoms.get_all_distances(mic=True) - get_contactMat(atoms, scale=1)
    count_contact = np.count_nonzero(diff <= max(tolerLims))/2
    count_tooClose = np.count_nonzero(min(tolerLims) >= diff)/2
    # # Count atoms pairs that are within a certain distance
    # tmpDiff = diff.copy()
    # tmpDiff[tmpDiff <= max(tolerLims)] = 1
    # tmpDiff[tmpDiff > max(tolerLims)] = 0

    # count_contact = tmpDiff.count_nonzero(min(tolerLims) <= tmpDiff <= max(tolerLims))
    # # Count atoms pairs that are too close to each other
    # tmpDiff = diff.copy()
    # tmpDiff[tmpDiff >= min(tolerLims)] = 0
    # tmpDiff[tmpDiff < min(tolerLims)] = 1
    # count_tooClose = tmpDiff.sum()
    # print(atoms, count_contact, count_tooClose)
    # print(diff)
    return count_contact - count_tooClose - len(diff)/2

# Main function for adding adsorbate fragments
def gen_solv_layer(
    interfc,
    solvList,
    xyzLims,
    toler=0.1
):
    tmpInterfc = interfc.copy()
    random.shuffle(solvList)
    ads_to_add = solvList.copy()
    n_added = 0
    n_trials = 0
    while n_added < len(ads_to_add):
        n_trials += 1
        tmpAds = ads_to_add[n_added].copy()
        del tmpAds[-1] # Here we delete the last atoms which must be dummy atom
        centreCoord = geom.rand_point_box(xyzLims)
        pos_ads = tmpAds.get_positions()
        # rotate along z
        rot = R.from_rotvec(- random.random()*np.pi * np.array([0, 0, 1]))
        pos_ads = rot.apply(pos_ads)
        # rotate randomly a bit
        rot_axis = geom.rand_direction()
        rot_angle = random.random()*np.pi/6
        rot = R.from_rotvec(- rot_angle * rot_axis)
        pos_ads = rot.apply(pos_ads)

        trj = []
        badStructure = True
        myContacts = count_bondContact(tmpAds, toler)
        myContacts_surf = count_bondContact(tmpInterfc.get_allAtoms(), toler)

        pos_ads_test = pos_ads + centreCoord
        tmpInterfc_test = tmpInterfc.copy()
        tmpAds_test = tmpAds.copy()
        tmpAds_test.set_positions(pos_ads_test)
        adsContacts_old = count_bondContact(
                tmpInterfc_test.get_adsAtoms(), toler)
        tmpInterfc_test.add_adsFrag(tmpAds_test)
        trj.append(tmpInterfc_test.get_allAtoms())
        surfContacts = count_bondContact(
                tmpInterfc_test.get_allAtoms(), toler)
        adsContacts = count_bondContact(
                tmpInterfc_test.get_adsAtoms(), toler)
        badContact = count_badContact(
                tmpInterfc_test.get_allAtoms())
        if badContact == count_badContact(interfc.get_allAtoms()) + count_badContact(tmpAds):
            if surfContacts == myContacts_surf+myContacts+1:
                if adsContacts == adsContacts_old + myContacts:
                    print(surfContacts, myContacts_surf, myContacts)
                    print(adsContacts, adsContacts_old, myContacts)
                    badStructure = False

        if not badStructure:
            tmpInterfc = tmpInterfc_test.copy()
            tmpInterfc.info = tmpInterfc_test.info
            n_added += 1
            print(f'{n_added}-{tmpAds.get_chemical_formula()}\tadded')
        write('test.xyz', trj)
        if n_trials >= 10000:
            print('FAIL! Exceeding 1000 attempts!')
            return None
    print('Acceptance rate: %.2f%%' % (len(solvList)/n_trials*100))
    return tmpInterfc
#        testInterfc.merge_adsorbate(Atoms(addElemList[ind_curr], [newAdsCoord]))

# Substrate
surf = Interface(
    './substrate.vasp',
    './substrate.vasp',
)
surf.print()

# Possible adsorbates
#methyl = read('./CH3.xyz')
co = read('./CO.xyz')
h = read('./H.xyz')
#h2o = read('./H2O.xyz')

# Add adsorbates to substrate
with connect('ini.db') as db:
    #db.write(surf.get_allAtoms())
    i = 0
    while i < 100:
        print(f'\n#{i}')
        mysurf = surf.copy()
        mysurf.rattle(0.5, zEnhance=True)
        adsList = random.randint(4,8)*[co] + random.randint(4,8)*[h]
        random.shuffle(adsList)
        slab_solvated = gen_solv_layer(
            mysurf, adsList,
            np.array([[0, 10.225], [0, 8.855], [12, 20]]), [-0.4,0.4])
   
        if slab_solvated is not None:
            slab_solvated.sortAds_frag()
            #slab_solvated.rattleMut()
            slab_solvated.preopt_lj()
            #slab_solvated.growMut_box(['CO'],np.array([[0, 15.336], [0, 15.336], [15, 20]]))
            #slab_solvated.growMut(['CO'])
            #slab_solvated.permuteMut()
            #slab_solvated.leachMut(['CO'])
            #slab_solvated.transMut()
            #slab_solvated.moveMut(['CO'])
            db.write(slab_solvated.get_allAtoms(), adsFrags="{}".format(slab_solvated.get_fragList()))
            i += 1

