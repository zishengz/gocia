
from gocia.ga.popGrandCanonPoly import PopulationGrandCanonicalPoly
#from gocia.ga.popGrandCanon import PopulationGrandCanonical
from gocia.utils.vasp import pos2pot, do_multiStep_opt
#from gocia.geom import del_freeMol
#from ase.io import read, write
import os
import numpy as np
import input

pop = PopulationGrandCanonicalPoly(
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
        xyzLims=np.array([[0,10.225],[0,8.85],[12,20]]),
        bondRejList = [['Cu','O']],#, ['O', 'O']],
        constrainTop=False,
        #moveOn=False,
        transVec=[[-2,2,-4,4],[-2,2,-4,4]]
        )
    # # Force odd coverage
    # if kid.get_adsAtoms().get_chemical_symbols().count('H')%2 == 0:
    #     kid = None
kid.preopt_hooke(cutoff=1.1, toler=0.1)
kid.print()
kid.write('POSCAR')

# VASP INPUT PREPARATION
pos2pot(input.pp_path)
os.system('cp POSCAR inp.vasp')
os.system('cp ../KPOINTS .')

# VASP CALCULATIONS
do_multiStep_opt(
    3, input.vasp_cmd,
    chkMol=True,
    list_keep=list(range(64)),
    zLim=input.zLim,
    potPath=input.pp_path,
    )

# ADD DATA& ClEAN UP
pop.add_vaspResult()
pop.natural_selection()
#os.system('rm WAVECAR CHG CHGCAR vasprun.xml POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')



