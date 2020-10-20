
import os, sys
import ase.io as fio
import ase.db as db
from gocia.interface import Interface
import gocia.utils.report as rp

# $ python plotCPK.py asym-4l.vasp Ga3N3.vasp
# # Timing: 1.4265 s, ~10% faster than BS
subsName = sys.argv[1]
surfName = sys.argv[2]
baseName = surfName.split('/')[-1].split('.')[0]

surf = Interface(
    tags = baseName,
    subAtoms = fio.read(subsName),
    allAtoms = fio.read(surfName),
    )
surf.info['eV'] = 0.1
surf.print()

surf.draw('CPK', title='%s, %.3f eV'%(surf.get_formula(), surf.info['eV']))
surf.draw('BS', title='%s, %.3f eV'%(surf.get_formula(),surf.info['eV']))
