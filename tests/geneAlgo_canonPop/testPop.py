
from gocia.ga.popCanon import PopulationCanonical
from ase.io import read, write
from ase.db import connect

pop = PopulationCanonical(gadb='rawInp.db', substrate='substrate.vasp', zLim=[7.5, 10.5])
pop.initializeDB()
aliveList = pop.get_ID('alive=1')
print(pop.get_valueOf('eV', aliveList))
print(pop.get_GMrow().eV)

with connect('gen1kid.db', append=False) as tmpDB:
    for i in range(20):
        kid = pop.gen_offspring()
        kid.preopt_hooke(cutoff=1.2, toler=0.1)
        tmpDB.write(kid.get_allAtoms())


