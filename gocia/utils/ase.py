


from ase.io import read, write
import os
import random
from gocia.geom import get_fragments, del_freeMol, is_bonded, detect_bond_between_adsFrag
from gocia.interface import Interface
from gocia.geom.frag import *

def geomopt_simple(atoms, my_calc, optimizer='LBFGS', fmax=0.05, label=None, fn_bkup=None, tmp_json='TMP.json'):
    atoms_opt = atoms.copy()
    atoms_opt.calc = my_calc
    write('ini.vasp', atoms)

    if label is not None:
        cwd = os.getcwd()
        try:
            os.mkdir(label)
        except:
            pass
        os.chdir(label)

    if optimizer is None:
        print('LOCAL OPTIMIZATION w/ built-in optimizer!')
        atoms_opt.get_potential_energy()
        if fn_bkup is not None:
            os.rename(f'{fn_bkup}', f'opt__{fn_bkup}')
    else:
        print(f'LOCAL OPTIMIZATION w/ {optimizer}')
        if optimizer == 'LBFGS':
            from ase.optimize import LBFGS
            dyn = LBFGS(atoms_opt, maxstep=0.1, trajectory=f'opt.traj', logfile=f'opt.log')

        if label is None:
            print(f'Optimizing {atoms_opt.get_chemical_formula()}')
        else:
            print(f'Optimizing {atoms_opt.get_chemical_formula()} in {label}')
            
        for i in range(1000):
            dyn.run(fmax=fmax, steps=1)
            if atoms_opt.get_forces().max() > 1000:
                print(f"Force exceeds threshold.")
                return None    

    if label is not None:
        os.chdir(cwd)

    # to solve some calculator issues
    if tmp_json:
        write(tmp_json, atoms_opt)
    atoms_opt = read(tmp_json)
    os.remove(tmp_json)

    return atoms_opt


def geomopt_multi(atoms, list_calc, optimizer='LBFGS', list_fmax=None, label=None, fn_bkup=None):

    atoms_opt = atoms.copy()
    if optimizer is not None and list_fmax is not None:
        if len(list_calc) != len(list_fmax):
            print('#calculators is inconsistent with #fmax!')
            exit()
        for i_step in range(len(list_calc)):
            print(f'STEP {i_step+1}')
            atoms_opt = geomopt_simple(atoms_opt, list_calc[i_step], optimizer, list_fmax[i_step], label=label)
            if label is None:
                label = '.'
            os.rename(f'{label}/opt.traj', f'{label}/opt{i_step+1}.traj')
    else:
        for i_step in range(len(list_calc)):
            print(f'STEP {i_step+1}')
            atoms_opt = geomopt_simple(atoms_opt, list_calc[i_step], optimizer, label=label)
            if fn_bkup is not None:
                if label is None:
                    label = '.'
                os.rename(f'{label}/{fn_bkup}', f'{label}/opt{i_step+1}__{fn_bkup}')
    return atoms_opt


def geomopt_iterate(atoms, my_calc, optimizer='LBFGS', fmax=None, label=None, chkMol=False, zLim=None, substrate='../substrate.vasp', fn_frag='fragments', list_keep=[0], has_fragList=False, has_fragSurfBond=False, check_rxn_frags=False, rmAtomsNotInBond=[], fn_bkup=None):
    # FRAGMENT-RELATED FUNCTIONS ARE NOT FINISHED YET
    if os.path.isfile(fn_frag):
        if read_frag(fn=fn_frag) is not None:
            has_fragList = True

    if type(substrate) is str:
        my_subs = read(substrate)
    else:
        my_subs = substrate
        
    continueRunning = True
    counter = 0
    
    opt_atoms = atoms.copy()
    while continueRunning:
        
        if opt_atoms is None:
            return None
        
        print(f'ITERATIVE OPTIMIZATION -- CYCLE: {counter + 1}')
        if optimizer is not None and fmax is not None:
            if type(my_calc) is list and type(fmax) is list:
                opt_atoms = geomopt_multi(opt_atoms, my_calc, optimizer, fmax, label=label, fn_bkup=fn_bkup)
            elif type(my_calc) is not list and type(fmax) is not list:
                opt_atoms = geomopt_simple(opt_atoms, my_calc, optimizer, fmax, label=label, fn_bkup=fn_bkup)
            else:
                print('Inconsistent #fmax and #calc\\ -- they should either be both lists ot both single objects/variables.')

            for fb in [f for f in os.listdir('.') if f[:3] == 'opt' and f[5:]=='.traj']:
                os.rename(fb, f'cyc{str(counter+1).zfill(2)}__{fb}')
        else:
            if type(my_calc) is list:
                opt_atoms = geomopt_multi(opt_atoms, my_calc, optimizer, label=label, fn_bkup=fn_bkup)
            elif type(my_calc) is not list:
                opt_atoms = geomopt_simple(opt_atoms, my_calc, optimizer, label=label, fn_bkup=fn_bkup)
            # back up files
            if fn_bkup is not None:
                if label is None:
                    label = '.'
                for fb in [f for f in os.listdir(label) if f[:3] == 'opt' and fn_bkup in f]:
                    os.rename(f'{label}/{fb}', f'{label}/cyc{str(counter+1).zfill(2)}__{fb}')
        
        if opt_atoms is None:
            return None
        
        natoms_old = len(opt_atoms)
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
            list1 = list(range(len(my_subs)))
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
            if len(geom_tmp) < len(opt_atoms):
                opt_atoms = geom_tmp
                if has_fragList and len(list_del) > 0:
                    update_frag_del(list_del, fn=fn_frag)

        # REMOVE ATOMS OUTSIDE THE SAMPLING BOX
        if zLim is not None:
            surf = Interface(
                opt_atoms,
                my_subs,
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

        natoms_new = len(opt_atoms)


        natoms_removed = natoms_old - natoms_new
        if natoms_removed > 0:
            print(f'Redo the last opt step due to removal of {natoms_removed} atoms')
            opt_atoms = opt_atoms.copy()
            continue
        elif natoms_removed == 0:
            print('The geometry is good! CONVERGED!')
            continueRunning = False

    # it still does not return a "good" calculator attribute
    # need to look more into this.
    if opt_atoms.calc.results is not None:
        # write(f'{label}_inner.json', opt_atoms)
        return opt_atoms

