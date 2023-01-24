
from ase.io import read, write
import os
from gocia.geom import get_fragments, del_freeMol
from gocia.interface import Interface
from gocia.geom.frag import *


def geomopt_DP(atoms, dp_file=None, fmax=0.1, optimizer='LBFGS'):
    from deepmd.calculator import DP
    calc_dp = DP(model=dp_file)
    atoms.calc = calc_dp
    write('POSCAR', atoms)
    if optimizer == 'LBFGS':
        from ase.optimize import LBFGS
        dyn = LBFGS(atoms, maxstep=0.05, trajectory='opt.traj')
    dyn.run(fmax=fmax)
    write('CONTCAR', atoms)
    return atoms

def do_multiStep_opt(atoms, dp_file=None, list_fmax=[0.1, 0.1], chkMol=False, zLim=None, substrate='../substrate.vasp', fn_frag='fragments', list_keep=[0], poscar='POSCAR'):
    step = len(list_fmax)
    write('POSCAR', atoms)
    continueRunning = True
    for i in range(step):
        if continueRunning:
            print(f'Optimization step: {i+1}')
            geomopt_DP(read('POSCAR'), dp_file=dp_file, fmax=list_fmax[i])
            os.system(f'cp CONTCAR out-{i}.vasp')
            if not chkMol:
                os.system('cp CONTCAR POSCAR')
            else:
                geom_tmp, list_del = del_freeMol(read('CONTCAR'), list_keep=list_keep)
                write('POSCAR', geom_tmp)
                my_fragList = read_frag(fn=fn_frag)
                if my_fragList is not None:
                    update_frag_del(list_del, fn=fn_frag)
                # Make sure the final structure has no free molecule
                if i == step:
                    if len(list_del) > 1:
                        os.system('touch BADSTRUCTURE')
                        continueRunning = False
                        continue
            if zLim is not None:
                surf = Interface(
                    read('POSCAR'),
                    substrate,
                    zLim = zLim
                )
                if surf.has_outsideBox():
                    # TODO: This breaks fragments
                    # need to modify
                    my_fragList = read_frag(fn=fn_frag)
                    if my_fragList is not None:
                        surf.del_outsideBox_frag()
                        write_frag(surf.get_fragList())
                    else:
                        surf.del_outsideBox()
                    surf.write('POSCAR')
                    if i == step:
                        os.system('touch BADSTRUCTURE')
                        continueRunning = False
                        continue
            my_fragList = read_frag(fn=fn_frag)
            if my_fragList is not None:
                struct = read('POSCAR')
                my_fragAtoms = [struct[f] for f in my_fragList]
                # Check connectivity
                list_del = []
                for i in range(len(my_fragList)):
                    if len(get_fragments(my_fragAtoms[i]))!=1:
                        list_del += my_fragList[i]
                if len(my_fragList) > 0:
                    print('Remove broken fragments containing:', list_del)
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    write('POSCAR', struct)

