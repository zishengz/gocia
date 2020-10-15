
import ase.io as ai
from gaia.interface import Interface
from gaia.utils import geom

surf = Interface(
    tags = 'test system',
    allAtoms = ai.read('test.xyz'),
    subAtoms = ai.read('CONTCAR')
)

surf2 = surf.copy()

surf.print_info()
surf2.print_info()

print(surf2.get_topLayer(0.5))
geom.chk_bondlength(surf2.get_allAtoms())
geom.get_bondpairs(surf2.get_allAtoms(), 0.85)
print(geom.get_neighbors(surf2.get_allAtoms(), 1)[0])
surf2.remove_adatom(1)
surf2.preopt_lj()
ai.write('out.vasp',surf2.get_allAtoms())

