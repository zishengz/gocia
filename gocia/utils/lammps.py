
from ase.io import read, write
import numpy as np

def atoms2lmpData(atoms):
    if type(atoms) is str:
        atoms = read(atoms, index='-1')
    write('inp.lammps-data', atoms, atom_style='charge')
    write('inp.vasp', atoms)

def get_last_frame(trajXyz, inp='inp.vasp'):
    s = read(inp)
    data = open(trajXyz).readlines()
    natoms = eval(data[3])
    pos_new = []
    for i in data[-natoms:]:
        pos_new.append([eval(j) for j in i.split()[1:]])
    s.set_positions(np.array(pos_new))
    return s

def get_ene(lmpOut='lmp.out'):
    data = open(lmpOut).readlines()
    for i in range(len(data)):
        if 'Energy initial, next-to-last, final =' in data[i]:
            data = data[i+1]
            break
    # return in eV
    return eval(data.split()[-1])/23.0605419
