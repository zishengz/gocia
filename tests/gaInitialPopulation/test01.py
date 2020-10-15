
import ase.io as ai
import ase.db as db
from gaia.interface import Interface
from gaia.utils import geom
from gaia.ga import build

surf = Interface(
    tags = 'test system',
    allAtoms = ai.read('CONTCAR'),
    subAtoms = ai.read('CONTCAR')
)

surf.print_info()

for i in range(50):
    newsurf = build.grow_adatom(
        surf,
        'N N'.split(),
        zLim=surf.zLim,
        ljopt=True
    )
    #newsurf.print_info()
    dtbs = db.connect('tmp.db')
    dtbs.write(newsurf.get_allAtoms(), {'eV':(250.213)})
    #print(dtbs)
    newsurf.write('out.json')
