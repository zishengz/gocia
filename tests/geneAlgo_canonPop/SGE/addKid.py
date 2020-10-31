from gocia.ga.popCanon import PopulationCanonical

pop = PopulationCanonical(gadb='../B4O6H1.db', substrate='../substrate.vasp', zLim=[7.5, 10.5])

pop.add_vaspResult()
pop.natural_selection()
