from gocia.ga.popCanon import PopulationCanonical

pop = PopulationCanonical(gadb='../B4O6H1.db', substrate='../substrate.vasp', zLim=[7.5, 10.5])

kid = None
while kid is None:
    kid = pop.gen_offspring()

kid.preopt_hooke(cutoff=1.1, toler=0.1)
kid.write('POSCAR')

