

import sys, os
import numpy as np
import ase.io as fio
from ase.db import connect
from gocia.utils.dbio import get_traj

def get_sorted_dist_list(atoms):
    """ Utility method used to calculate the sorted distance list
        describing the cluster in atoms. """
    numbers = atoms.numbers
    unique_types = set(numbers)
    pair_cor = dict()
    for n in unique_types:
        i_un = [i for i in range(len(atoms)) if atoms[i].number == n]
        d = []
        for i, n1 in enumerate(i_un):
            for n2 in i_un[i + 1:]:
                d.append(atoms.get_distance(n1, n2, mic=True))
        d.sort()
        pair_cor[n] = np.array(d)
    return pair_cor

def srtDist_similar(a1, a2, delta_rel=0.001, d_max=3):
    """ Private method for calculating the structural difference. """
    p1 = get_sorted_dist_list(a1)
    p2 = get_sorted_dist_list(a2)
    numbers = a1.numbers
    total_cum_diff = 0.
    max_diff = 0
    for n in p1.keys():
        cum_diff = 0.
        c1 = p1[n]
        c2 = p2[n]
        assert len(c1) == len(c2)
        if len(c1) == 0:
            continue
        t_size = np.sum(c1)
        d = np.abs(c1 - c2)
        cum_diff = np.sum(d)
        max_diff = np.max(d)
        ntype = float(sum([i == n for i in numbers]))
        total_cum_diff += cum_diff / t_size * ntype / float(len(numbers))
    print(total_cum_diff, cum_diff)
    if total_cum_diff < delta_rel and max_diff < d_max:
        return True
    else:
        return False

def comp_srtDist(a1, traj, delta_rel=0.001, d_max=3):
    simList = [0]*len(traj)
    for i in range(len(traj)):
        if srtDist_similar(a1, traj[i], delta_rel=delta_rel, d_max=d_max):
            simList[i] = 1
            print('it looks like %i'%(i+1))
    return simList

def srtDist_similar_zz(a1, a2, delta_rel=1e-4, d_max=0.1):
    p1 = a1.get_all_distances(mic=True).flatten()
    p2 = a2.get_all_distances(mic=True).flatten()
    p1, p2 = np.sort(p1), np.sort(p2)
    cum_diff = np.abs(p1 - p2)
    total_cum_diff = cum_diff.sum() / 2 / (p1 + p2).sum()
    max_diff = cum_diff.max()
    print(total_cum_diff, max_diff)
    if total_cum_diff < delta_rel and max_diff < d_max:
        return True
    else:
        return False

def comp_srtDist_zz(a1, traj, delta_rel=1e-4, d_max=0.1):
    simList = [0]*len(traj)
    for i in range(len(traj)):
        if srtDist_similar_zz(a1, traj[i], delta_rel=delta_rel, d_max=d_max):
            simList[i] = 1
            print('it looks like %i'%(i+1))
    return simList

def eigen_similar_zz(a1, a2, delta_rel=0.01, d_max=1):
    p1 = a1.get_all_distances(mic=True)
    p2 = a2.get_all_distances(mic=True)
    p1 = 1/np.linalg.eigh(p1)[0]
    p2 = 1/np.linalg.eigh(p2)[0]
    cum_diff = np.abs(p1 - p2)
    total_cum_diff = abs( cum_diff.sum() / 2 / (p1 + p2).sum() )
    max_diff = cum_diff.max()
    print(total_cum_diff, max_diff)
    if total_cum_diff < delta_rel and max_diff < d_max:
        return True
    else:
        return False

def eigen_srtDist_zz(a1, traj, delta_rel=0.01, d_max=1):
    simList = [0]*len(traj)
    for i in range(len(traj)):
        if eigen_similar_zz(a1, traj[i], delta_rel=delta_rel, d_max=d_max):
            simList[i] = 1
            print('it looks like %i'%(i))
    return simList

trajDBName = sys.argv[1]
geomCutoff = 0.1
enerCutoff = 0.1

print(' * Extracting data from database: %s'%trajDBName)
rawDB = connect(trajDBName)
traj = get_traj(rawDB.select())

comp_srtDist_zz(traj[0], traj[:30])


