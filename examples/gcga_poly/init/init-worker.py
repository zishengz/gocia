from gocia.utils.vasp import pos2pot, do_multiStep_opt
from gocia.geom import get_fragments
from ase.io import read
import os
import numpy as np
import input

# VASP INPUT PREPARATION
pos2pot(input.pp_path)
os.system('cp POSCAR inp.vasp')
os.system('cp ../KPOINTS .')

# 3-step VASP CALCULATIONS
do_multiStep_opt(
    3, input.vasp_cmd,
    chkMol=True,
#    list_keep=list(range(64)),
#    zLim=input.zLim,
    potPath=input.pp_path,
    )

# Simple 1-step GEOM OPT
#os.system('cp ../INCAR .')
#os.system(input.vasp_cmd)

