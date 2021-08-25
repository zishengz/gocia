from gocia.utils.vasp import pos2pot
from gocia.geom import get_fragments
from ase.io import read
import os
import numpy as np
import input

# VASP INPUT PREPARATION
pos2pot(input.pp_path)
os.system('cp POSCAR inp.vasp')
os.system('cp ../KPOINTS .')

# VASP CALCULATIONS
for i in range(1,3):
    os.system('../INCAR-%i INCAR'%i)
    os.system(input.vasp_cmd)
    os.system('cp CONTCAR out-%i.vasp'%i)
    if len(get_fragments(read('CONTCAR'))) > 1:
        os.system('touch BADSTRUCTURE')
        exit()

os.system('rm WAVECAR CHG CHGCAR vasprun.xml POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')



