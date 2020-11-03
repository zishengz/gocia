from gocia.interface import Interface
import sys
from gocia.ga.crossover import crossover_snsSurf_2d_GC

subs = sys.argv[1]
surf1 = sys.argv[2]
surf2 = sys.argv[3]

#surf0 = Interface(subAtoms=subs, allAtoms=subs)
surf1 = Interface(allAtoms=surf1, subAtoms=subs)
surf2 = Interface(allAtoms=surf2, subAtoms=subs)

surf1.print()
surf2.print()

#print(surf1.get_fixBufPos())
#print(surf1.get_fixBufAtoms())

nova = crossover_snsSurf_2d_GC(surf1, surf2)
nova.print()
nova.write('nova.vasp')
