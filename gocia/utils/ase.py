

from ase.io import read, write
import os
import random
from gocia.geom import get_fragments, del_freeMol, is_bonded, detect_bond_between_adsFrag
from gocia.interface import Interface
from gocia.geom.frag import *

def geomopt_simple(atoms, my_calc, fmax=0.1, label=None, optimizer='LBFGS'):
    atoms.calc = my_calc

    if label is not None:
        cwd = os.getcwd()
        try:
            os.mkdir(label)
        except:
            pass
        os.chdir(label)

    if optimizer == 'LBFGS':
        from ase.optimize import LBFGS
        dyn = LBFGS(atoms, maxstep=0.05, trajectory=f'opt.traj', logfile=f'opt.log')

    print(f'Optimizing {atoms.get_chemical_formula()} in {label}')
    dyn.run(fmax=fmax)

    if label is not None:
        os.chdir(cwd)
    return atoms


def geomopt_iterate(atoms, my_calc, fmax=0.1, label='test', optimizer='LBFGS', chkMol=False, zLim=None, substrate='../substrate.vasp', fn_frag='fragments', list_keep=[0], has_fragList=False, has_fragSurfBond=False, check_rxn_frags=False, rmAtomsNotInBond=[]):
    # FRAGMENT-RELATED FUNCTIONS ARE NOT FINISHED YET
    if fn_frag in os.listdir():
        if read_frag(fn=fn_frag) is not None:
            has_fragList = True
        
    continueRunning = True
    counter = 0
    my_atoms = atoms.copy()
    while continueRunning:
            print(f'Optimization cycle: {counter}')
            opt_atoms = geomopt_simple(my_atoms, my_calc, fmax=fmax, label=label, optimizer=optimizer)
            counter += 1

            # REMOVE BROKEN FRAGMENTS
            if has_fragList:
                my_fragList = read_frag(fn=fn_frag)
                struct = opt_atoms.copy()
                my_fragAtoms = [struct[f] for f in my_fragList]
                # Check connectivity
                list_del = []
                for i in range(len(my_fragList)):
                    if len(get_fragments(my_fragAtoms[i]))!=1:
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                if len(list_del) > 0:
                    print('Remove broken fragments containing:', list_del)
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    opt_atoms = struct

            # REMOVE FRAGMENTS THAT DOES NOT BIND TO SURFACE ATOMS
            # E.G. [SURF]-[C-O]-[C-O]
            #      the latter CO is bonded to an adorbate but not the surface
            if has_fragList and has_fragSurfBond:
                my_fragList = read_frag(fn=fn_frag)
                struct = opt_atoms.copy()
                list1 = list(range(len(read(substrate))))
                list_del = []
                for i in range(len(my_fragList)):
                    if not is_bonded(struct, list1, my_fragList[i]):
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                if len(list_del) > 0:
                    print('Remove fragments not bonded to the surface:', list_del)
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    opt_atoms = struct

            # REMOVE FRAGMENTS THAT FROM BOND WITH OTHER FRAGMENTS
            if has_fragList and check_rxn_frags:
                my_fragList = read_frag(fn=fn_frag)
                struct = opt_atoms.copy()
                ads_bonds = detect_bond_between_adsFrag(struct, my_fragList)
                if len(ads_bonds) > 0:
                    list_del = []
                    for b in ads_bonds:
                        i = random.choice(b[0:2])
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                    if len(list_del) > 0:
                        print('Remove associated fragments containing:', list_del)
                        update_frag_del(list_del, fn=fn_frag)
                        del struct[list_del]
                        opt_atoms = struct

            # REMOVE FREE MOLECULES DESORBED FROM THE SURFACE
            if chkMol:
                geom_tmp, list_del = del_freeMol(opt_atoms, list_keep=list_keep)
                opt_atoms = geom_tmp
                if has_fragList and len(list_del) > 0:
                    update_frag_del(list_del, fn=fn_frag)

            # REMOVE ATOMS OUTSIDE THE SAMPLING BOX
            if zLim is not None:
                surf = Interface(
                    opt_atoms,
                    substrate,
                    zLim = zLim
                )
                if surf.has_outsideBox():
                    if has_fragList:
                        surf.del_outsideBox_frag(fn_frag)
                    else:
                        surf.del_outsideBox()
                    opt_atoms = surf.get_allAtoms()

            # REMOVE ATOMS THAT DO NOT FORM SPECIFIED BONDS
            if len(rmAtomsNotInBond) > 0:
                # e.g. if [['H', 'Pt']]
                # 'H' will be removed if not in Pt-H
                struct = opt_atoms.copy()
                list_del = []
                for p in rmAtomsNotInBond:
                    list_a0 = [a.index for a in struct if a.symbol==p[0]]
                    list_a1 = [a.index for a in struct if a.symbol==p[1]]
                    for a0 in list_a0:
                        if not is_bonded(struct, [a0], list_a1):
                            list_del.append(a0)
                if len(list_del) > 0:
                    print(f'Atoms {list_del} are removed becasue not in required bonds.')
                    del struct[list_del]
                    opt_atoms = struct

            # REMOVE BROKEN FRAGMENTS
            if has_fragList:
                my_fragList = read_frag(fn=fn_frag)
                struct = opt_atoms.copy()
                my_fragAtoms = [struct[f] for f in my_fragList]
                # Check connectivity
                list_del = []
                for i in range(len(my_fragList)):
                    if len(get_fragments(my_fragAtoms[i]))!=1:
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                if len(list_del) > 0:
                    print('Remove broken fragments containing:', list_del)
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    opt_atoms = struct

            natoms_removed = len(my_atoms) - len(opt_atoms)
            my_atoms = opt_atoms
            if natoms_removed > 0:
                print(f'Redo the last opt step due to removal of {natoms_removed} atoms')
                continue
            elif natoms_removed == 0:
                continueRunning = False
    return my_atoms
