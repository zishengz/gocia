from ase.db import connect
from gocia.ga.crossover import crossover_snsSurf_2d_GC_poly
from gocia.interface import Interface
from gocia.ga.popGrandCanonPoly import PopulationGrandCanonicalPoly

###
#   practiceDB.py
#   
#   This script is to practice working with the database of structures, 
#   specifically getting the list of atom indices in each adsorbate fragment
#   and using this info to generate kid structures from parents in the database
#
#   @author: Winston Gee
# 
###


db = connect('ini.db')

surf = Interface(
    './substrate.vasp',
    './substrate.vasp'
)

for i in list(range(2,5)):
    row = db.get(id=i)
    print("{}".format(row.key_value_pairs))
    print("{0} : {1} : {2}".format(row.id,row.formula,row.get('adsFrags')))

    pop = PopulationGrandCanonicalPoly(
        gadb='ini.db',
        substrate='substrate.vasp',
        popSize = 5,
        convergeCrit=50,
        subsPot = -885.44371, # not true but just a value to try for now
        chemPotDict = {'CO':-14.4},
        zLim = [15, 20]
        )

    a1 = db.get(id=i-1).toatoms()
    a2 = db.get(id=i).toatoms()
    f1 = db.get(id=i-1).get('adsFrags')
    f2 = db.get(id=i).get('adsFrags')
    info1 = pop.convertFragStrToInfo(f1)
    info2 = pop.convertFragStrToInfo(f2)
    print(f1, info1, f2, info2)
    surf1 = Interface(a1, surf.get_allAtoms(), info=info1)
    surf2 = Interface(a2, surf.get_allAtoms(), info=info2)

    kid = crossover_snsSurf_2d_GC_poly(surf1,surf2,tolerance=0.75)
    print(kid.get_fragList())
    kid.print()
    #kid.moveMut_frag(['CO'])

    db.write(kid.get_allAtoms(), adsFrags="{}".format(kid.get_fragList()))
