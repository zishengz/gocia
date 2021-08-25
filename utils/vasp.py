import os

def get_elems(poscar='POSCAR', potDict=None):
    elems = open(poscar).readlines()[5].split()
    if potDict is not None:
        elems = [potDict[e] for e in elems]
    return elems

def gen_POTCAR(potPath, elems):
    with open('POTCAR', 'w') as f:
        for e in elems:
            f.write(open('%s/%s/POTCAR'%(potPath, e)).read())

def pos2pot(potPath, poscar='POSCAR', potDict=None):
    gen_POTCAR(potPath, get_elems(poscar, potDict))

