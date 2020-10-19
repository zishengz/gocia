
import numpy as np

def get_keyVal(dbobj, keyVal):
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