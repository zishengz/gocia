
import ase.io as fio
from gocia.interface import Interface

from ase.calculators.emt import EMT
from ase.optimize.bfgs import BFGS

subs = fio.read('subs-ini.vasp')
subs.calc = EMT()
geomOpt = BFGS(subs, maxstep=0.05)
geomOpt.run(fmax=0.001, steps=200)
fio.write('subs-opt.vasp', subs)

