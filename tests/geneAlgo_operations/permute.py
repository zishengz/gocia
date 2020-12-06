import sys
sys.path.append('/home/zisheng/windir/github')
from ase.io import read, write
import numpy as np
from gocia.geom import cart2frac, frac2cart
from gocia.interface import Interface

surf = Interface(
    'test.vasp',
    'substrate.vasp'
)
surf.print()

surf.permuteMut()

surf.print()
surf.write('permute.vasp')





