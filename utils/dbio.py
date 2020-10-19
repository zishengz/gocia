
import numpy as np

def get_keyVal(dbRows, keyVal):
    '''
    Returns 1-D list
    '''
    tmp = []
    for r in dbRows:
        tmp.append(r[keyVal])
    return tmp

def get_traj(dbRows):
    tmp = []
    for r in dbRows:
        a = r.toatoms()
        a.info = r.key_value_pairs
        tmp.append(a)
    return tmp
