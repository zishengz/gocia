import os
import numpy as np
from ase.io import read, write
from ase.db import connect
from gocia.utils.ase import geomopt_iterate
from gocia.geom.frag import read_frag
from xtb.ase.calculator import XTB
from gocia.utils.dbio import get_traj, vasp2db, get_projName
from gocia.ensemble import clusterIsomerAbs
from ase.constraints import FixAtoms
import input
from ase.neighborlist import NeighborList, natural_cutoffs

inp = 'ini.db'
structs_raw = read('ini.db', ':')
structs_sub = read('sub-POSCAR')

structs_opt = []

mycalc = XTB(method="GFN2-xTB")

with connect('gcga.db') as opt_db:
    opt_db.write(
    s_new,
    eV = s_new.get_potential_energy(),
    adsFrags='%s' % (fn_frag),
    done=1,
    )

accept_list = []
reject_list = []
with connect(inp) as db:
    for i in range(len(db)):
        my_label = f's{str(i).zfill(6)}'
    
        # OPTIMIZE THE STRUCTURE
        # Write fragment to file
        frag_filename = f'{my_label}.fragments'
        with open(frag_filename, 'w') as f:
            f.write(db.get(id=i+1).get('adsFrags'))

        fmax = 0.02
        s_new = geomopt_iterate(
                structs_raw[i],
                my_calc, fmax=fmax,
                label=my_label,
                relax_steps=500,
                chkMol=True,
                zLim=input.zLim,
                substrate='sub-POSCAR',
                list_keep=list(range(64)),
                fn_frag='%s.fragments' % (my_label),
                has_fragSurfBond=True,
                check_rxn_frags=True,
                )
        s_new.calc = my_calc
        max_force = np.max(np.abs(s_new.get_forces().flatten()))

        fn_frag = read_frag(fn=frag_filename)
        print(my_label, fn_frag)
        
        if not examine_O2_molecule_presents(s_new) and not examine_CO_molecule_presents(s_new) and len(fn_frag) != 0 and max_force <= fmax:
            accept_list.append(i)
            with connect('gcga.db') as opt_db:
                opt_db.write(
                    s_new,
                    eV = s_new.get_potential_energy(),
                    adsFrags='%s' % (fn_frag),
                    done=1,
                )
        else:
            reject_list.append(i)

with connect('gcga.db') as rawDB:
    traj = get_traj(rawDB.select())
    clusterIsomerAbs(
        traj,
        eneToler=0.05,
        geomToler1=1e-3,
        geomToler2=0.5,
        outName='sorted'
    )

os.makedirs("z_init_pop")
os.system("mv s0* z_init_pop/")
os.system("cp sorted.db z_init_pop/")
os.system("cp gcga.db z_init_pop/")
print(f"rejected structures={reject_list}")
print(f"accepted structures={accept_list}")
