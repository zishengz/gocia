from ase.io import read, write
from ase.optimize import FIRE, LBFGS
from ase.db import connect
from mace.calculators import MACECalculator
import sys
import os
from ase.db import connect
from ase import Atoms # Import Atoms
import numpy as np # Import numpy
from tqdm import tqdm

from gocia.interface import Interface

# Check model path. MACE foundation model as an example
path_model = 'xxxxxx/mace-mpa-0-medium.model'
if not os.path.exists(path_model):
    raise FileNotFoundError(f"Model not found: {path_model}")


def init_calc():
    """Initializes and returns a MACE calculator instance."""
    return MACECalculator(model_paths=path_model, device='cpu', default_dtype='float32', dtype='float32')#, E0s='none')

def optimize_structure(atoms, init_calc, fmax=0.2, optimizer=LBFGS, logfile=None):
    """Optimizes the geometry of a single ASE Atoms object."""
    atoms.calc = init_calc()
    optimizer = optimizer(atoms, maxstep=0.05)
    for i in range(1000):
        optimizer.run(fmax=fmax, steps=1)
        if atoms.get_forces().max() > 500:
            print(f"Force exceeds threshold.")
            return None     
        elif optimizer.converged():
            break
    return atoms

# Read input structure
inpname = sys.argv[1]
basename = os.path.splitext(inpname)[0]

counter_good = 0


with connect(inpname) as mydb:
    for i, r in enumerate(mydb.select()):
        
        if f's{str(i).zfill(6)}.json' in os.listdir('.'):
            continue
        
        print(i, r.toatoms().get_chemical_formula())
        
        if r.get('optimized', False):
            continue
        
        

        struct = r.toatoms()
        
        # surf = Interface(struct, struct)
        # surf.preopt_hooke(cutoff = 1.0, toler = 0.05)
        
        struct = optimize_structure(struct, init_calc)

        if struct is not None:
            
            write(f's{str(i).zfill(6)}.json', struct)

            # Update the database entry, explicitly passing custom keys for energy and forces
            mydb.update(id = r.id, atoms = read(f's{str(i).zfill(6)}.json'), optimized = True)

            os.remove(f's{str(i).zfill(6)}.json')
            
            counter_good += 1

            # os.remove('tmp.json')
        
