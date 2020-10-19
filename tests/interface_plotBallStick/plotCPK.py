
import os, sys
import ase.io as fio
import ase.db as db
from gocia.interface import Interface
from gocia.utils import visualize

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
visualize.drawCPK(
    surf,
    outName=baseName
    )
