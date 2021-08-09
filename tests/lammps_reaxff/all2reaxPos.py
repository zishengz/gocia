from ase.io import read, write
import sys

fn_inp = sys.argv[1]

s = read(fn_inp)

write('%s.lammps-data'%(fn_inp.split('.')[0]), s, atom_style='charge')
