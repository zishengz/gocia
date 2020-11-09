from gocia.interface import Interface
import sys, os
from gocia.ga.crossover import crossover_snsSurf_2d_GC

subs = sys.argv[1]
surf1 = sys.argv[2]
surf2 = sys.argv[3]

#surf0 = Interface(subAtoms=subs, allAtoms=subs)
surf1 = Interface(allAtoms=surf1, subAtoms=subs)
surf2 = Interface(allAtoms=surf2, subAtoms=subs)

surf1.print()
surf1.draw('CPK')
os.system('mv Ga47N51-cpk.png mom.png')

surf2.print()
surf2.draw('CPK')
os.system('mv Ga47N51-cpk.png dad.png')

#print(surf1.get_fixBufPos())
#print(surf1.get_fixBufAtoms())

nova = crossover_snsSurf_2d_GC(surf1, surf2)
nova.print()
nova.write('nova.vasp')
nova.draw('CPK')
os.system('mv Ga47N51-cpk.png kid.png')
tmp = nova.copy()

nova = tmp.copy()
nova.rattleMut()
nova.draw('CPK')
os.system('mv Ga47N51-cpk.png kid-rattle.png')

nova = tmp.copy()
nova.growMut(['N'])
nova.draw('CPK')
os.system('mv Ga47N51-cpk.png kid-grow.png')

nova = tmp.copy()
nova.leachMut(['N'])
nova.draw('CPK')
os.system('mv Ga47N51-cpk.png kid-leach.png')
