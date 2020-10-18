
import os, sys
import ase.io as fio
import ase.db as db
from gocia.interface import Interface
from gocia.utils import visualize

surfName = sys.argv[1]
baseName = surfName.split('/')[-1].split('.')[0]

# # With Interface object constructed
# # Timing: 1.6935 s
# surf = Interface(
#     tags = baseName,
#     allAtoms = fio.read(surfName),
#     subAtoms = fio.read(surfName)
#     )
# visualize.drawBSsurf(
#     surf.get_allAtoms(),
#     zBuffer=[13.1, 15.8],
#     pseudoBond=True,
#     outName=baseName
#     )

# Directly use Atoms object
# Timing: 1.6355 s (3.55% speed-up)
surf = fio.read(surfName)
visualize.drawBSsurf(
    surf,
    zBuffer=[13.1, 15.8],
    pseudoBond=True,
    outName=baseName
    )
