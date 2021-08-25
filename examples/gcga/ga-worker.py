from gocia.ga.popGrandCanon import PopulationGrandCanonical
from gocia.utils.vasp import pos2pot
from gocia.geom import get_fragments
from ase.io import read
import os
import numpy as np
import input

pop = PopulationGrandCanonical(
    gadb='../gcga.db',
    substrate='../substrate.vasp',
    popSize = input.popSize,
    subsPot = input.subsPot,
    chemPotDict = input.chemPotDict,
    zLim = input.zLim
    )

# Child Generation
kid = None
while kid is None:
    kid = pop.gen_offspring_box(
        mutRate=0.4,
        xyzLims=np.array([[0,11.25],[0,9.74],[8,13.0]]),
        bondRejList = [['H','H'], ['Ni','Ni'], ['Ni', 'H'], ['O', 'O']],
        constrainTop=True,
        transVec=[[-3,3],[-3,3]]
        )
kid.preopt_hooke(cutoff=1.2, toler=0.1)
kid.write('POSCAR')

# VASP INPUT PREPARATION
pos2pot(input.pp_path)
os.system('cp POSCAR inp.vasp')
os.system('cp ../KPOINTS .')

# VASP CALCULATIONS
for i in range(1,3):
    os.system('cp ../INCAR-%i INCAR'%i)
    os.system(input.vasp_cmd)
    os.system('cp CONTCAR out-%i.vasp'%i)
    if len(get_fragments(read('CONTCAR'))) > 1:
        os.system('touch BADSTRUCTURE')
        exit()
    else:
        os.system('cp CONTCAR POSCAR')

# ADD DATA& ClEAN UP
pop.add_vaspResult()
pop.natural_selection()
os.system('rm WAVECAR CHG CHGCAR vasprun.xml POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')



