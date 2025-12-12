import random
import numpy as np
from ase.units import kB


def pass_Metropilis(delta_e, temperature):
    if delta_e < 0:
        return True
    elif random.random() < np.exp(-delta_e / (kB * temperature)):
        return True
    else:
        False


