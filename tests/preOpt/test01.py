
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

newsurf = build.grow_adatom(
        surf,
        'N N N'.split(),
        toler=0.75,
        rattle=True, rattleStdev=0.25,
        zLim=surf.zLim
    )
newsurf.write('ini.xyz')

# surfOpt = newsurf.copy()
# surfOpt.preopt_lj(toler=0.1)
# surfOpt.write('lj.xyz')

# surfOpt = newsurf.copy()
# surfOpt.preopt_hooke(toler=0.1)
# surfOpt.write('hk.xyz')

for i in range(20):
    newsurf = build.grow_adatom(
        surf,
        'N N N'.split(),
        toler=0.75,
        rattle=True,
        zLim=surf.zLim
    )
    newsurf.preopt_hooke(toler = 0.1)
    newsurf.write('N3g0.db')

