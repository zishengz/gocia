import os
import numpy as np
import random
from ase.io import read, write
from gocia.geom import get_fragments, del_freeMol, is_bonded, detect_bond_between_adsFrag
from gocia.interface import Interface
from gocia.geom.frag import *


def get_elems(poscar='POSCAR', potDict=None):
    elems = open(poscar).readlines()[5].split()
    if potDict is not None:
        elems = [potDict[e] for e in elems]
    return elems


def gen_POTCAR(potPath, elems):
    with open('POTCAR', 'w') as f:
        for e in elems:
            f.write(open('%s/%s/POTCAR' % (potPath, e)).read())


def pos2pot(potPath, poscar='POSCAR', potDict=None):
    gen_POTCAR(potPath, get_elems(poscar, potDict))


def add_hubbard(elems, incar='INCAR', UDict=None):
    ldaul = [str(UDict[elem]['L']) for elem in elems]
    ldauu = [f"{UDict[elem]['U']:.3f}" for elem in elems]
    ldauj = [f"{UDict[elem]['J']:.3f}" for elem in elems]

    with open(incar, 'w') as f:
        f.write("LDAU = .TRUE.\n")
        f.write("LDAUL = " + " ".join(ldaul) + "\n")
        f.write("LDAUU = " + " ".join(ldauu) + "\n")
        f.write("LDAUJ = " + " ".join(ldauj) + "\n")
        f.write("LDAUPRINT = 0\n")
        f.write("LDAUTYPE = 2\n")
        f.write("LMAXMIX = 4\n")


def is_vaspSuccess(jobdir='.'):
    return 'E0' in open(f'{jobdir}/OSZICAR').readlines()[-1]


# TODO: DECOUPLE THE GEOMETRY CHECKER & MULTI_STEP OPT?

def do_multiStep_opt(step=3, vasp_cmd='', chkMol=False, zLim=None, substrate='../substrate.vasp', fn_frag='fragments', list_keep=[0], potPath=None, poscar='POSCAR', potDict=None, has_fragList=False, has_fragSurfBond=False, check_rxn_frags=False, rmAtomsNotInBond=[], addU=False, UDict=None,):
    if read_frag(fn=fn_frag) is not None:
        has_fragList = True
    continueRunning = True
    counter = 1
    while counter <= step:
        if continueRunning:
            print(f'Optimization step: {counter}')
            os.system(f'cp ../INCAR-{counter} INCAR')
            if addU:
                add_hubbard(get_elems(poscar), incar='INCAR', UDict=UDict)
            pos2pot(potPath, poscar=poscar, potDict=potDict)
            os.system(vasp_cmd)
            if not is_vaspSuccess():
                os.system('touch FAIL')
                exit()
            os.system(f'cp CONTCAR out-{counter}.vasp')
            os.system(f'cp vasprun.xml vasprun-{counter}.xml')
            os.system('cp CONTCAR POSCAR')
            atom_tmp = read('CONTCAR')
            counter += 1

            # REMOVE ATOMS OUTSIDE THE SAMPLING BOX
            if zLim is not None:
                surf = Interface(
                    read('POSCAR'),
                    substrate,
                    zLim = zLim
                )
                if surf.has_outsideBox():
                    if has_fragList:
                        surf.del_outsideBox_frag(fn_frag)
                    else:
                        surf.del_outsideBox()
                    surf.write('POSCAR')

            # REMOVE BROKEN FRAGMENTS
            if has_fragList:
                my_fragList = read_frag(fn=fn_frag)
                struct = read('POSCAR')
                my_fragAtoms = [struct[f] for f in my_fragList]
                # Check connectivity
                list_del = []
                for i in range(len(my_fragList)):
                    if len(get_fragments(my_fragAtoms[i]))!=1:
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                if len(list_del) > 0:
                    list_del.sort()
                    print('Remove broken fragments containing:', list_del)
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    write('POSCAR', struct)

            # REMOVE FRAGMENTS THAT DOES NOT BIND TO SURFACE ATOMS
            # E.G. [SURF]-[C-O]-[C-O]
            #      the latter CO is bonded to an adorbate but not the surface
            if has_fragList and has_fragSurfBond:
                my_fragList = read_frag(fn=fn_frag)
                struct = read('POSCAR')
                list1 = list(range(len(read(substrate))))
                list_del = []
                for i in range(len(my_fragList)):
                    if not is_bonded(struct, list1, my_fragList[i]):
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                if len(list_del) > 0:
                    list_del.sort()
                    print('Remove associated fragments containing:', list_del)
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    write('POSCAR', struct)

            # REMOVE FRAGMENTS THAT FROM BOND WITH OTHER FRAGMENTS
            if has_fragList and check_rxn_frags:
                my_fragList = read_frag(fn=fn_frag)
                struct = read('POSCAR')
                ads_bonds = detect_bond_between_adsFrag(struct, my_fragList)
                if len(ads_bonds) > 0:
                    list_del = []
                    for b in ads_bonds:
                        i = random.choice(b[0:2])
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                    if len(list_del) > 0:
                        list_del.sort()
                        print('Remove associated fragments containing:', list_del)
                        update_frag_del(list_del, fn=fn_frag)
                        del struct[list_del]
                        write('POSCAR', struct)

            # REMOVE FREE MOLECULES DESORBED FROM THE SURFACE
            if chkMol:
                geom_tmp, list_del = del_freeMol(read('POSCAR'), list_keep=list_keep)
                write('POSCAR', geom_tmp)
                if has_fragList and len(list_del) > 0:
                    list_del.sort()
                    update_frag_del(list_del, fn=fn_frag)

            # REMOVE ATOMS THAT DO NOT FORM SPECIFIED BONDS
            if len(rmAtomsNotInBond) > 0:
                # e.g. if [['H', 'Pt']]
                # 'H' will be removed if not in Pt-H
                struct = read('POSCAR')
                list_del = []
                for p in rmAtomsNotInBond:
                    list_a0 = [a.index for a in struct if a.symbol==p[0]]
                    list_a1 = [a.index for a in struct if a.symbol==p[1]]
                    for a0 in list_a0:
                        if not is_bonded(struct, [a0], list_a1):
                            list_del.append(a0)
                if len(list_del) > 0:
                    list_del.sort()
                    print(f'Atoms {list_del} are removed becasue not in required bonds.')
                    write('POSCAR', struct)
                    del struct

            # REMOVE BROKEN FRAGMENTS
            if has_fragList:
                my_fragList = read_frag(fn=fn_frag)
                struct = read('POSCAR')
                my_fragAtoms = [struct[f] for f in my_fragList]
                # Check connectivity
                list_del = []
                for i in range(len(my_fragList)):
                    if len(get_fragments(my_fragAtoms[i]))!=1:
                        for i_del in my_fragList[i]:
                            if i_del not in list_del:
                                list_del.append(i_del)
                if len(list_del) > 0:
                    list_del.sort()
                    print('Remove broken fragments containing:', list_del)
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    write('POSCAR', struct)

            natoms_removed = len(atom_tmp) - len(read('POSCAR'))
            if natoms_removed > 0 and counter > step:
                print(f'Redo the last opt step due to removal of {natoms_removed} atoms')
                counter -= 1 # redo the last opt step if something is removed
                continue
    os.system('rm WAVECAR CHG CHGCAR POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')


# Below are for surface charging calculations

def get_neu_nelect(poscar='POSCAR', potcar='POTCAR', shift=0):
    elem_list = [l.split()[-2]
                 for l in open(potcar).readlines() if 'TITEL' in l]
    elem_list = [e.split('_')[0] if '_' in e else e for e in elem_list]
    zval_list = [eval(l.split()[-4])
                 for l in open(potcar).readlines() if 'ZVAL' in l]
    atom_list = read(poscar).get_chemical_symbols()
    return sum([zval_list[elem_list.index(a)] for a in atom_list]) + shift


def extractVASPsol(dirName='.', shift=0):
    tmpHome = os.getcwd()
    os.chdir(dirName)
    val = get_neu_nelect(shift=shift)
    fermi_shift = eval([l for l in open('out').readlines()
                        if 'FERMI_SHIFT' in l][-1].split()[2])
    fermi_ene = eval([l for l in open('OUTCAR').readlines()
                      if 'E-fermi' in l][-1].split()[2])
    nelect = eval([l for l in open('INCAR').readlines()
                   if 'NELECT' in l][-1].split('=')[-1])
    energy = eval([l for l in open('OSZICAR').readlines() if 'F=' in l][-1].split()[4])
    os.chdir(tmpHome)
    return val, nelect, energy, fermi_ene, fermi_shift


def pb_calc(val, nelect, energy, fermi_ene, fermi_shift, u_ref=4.44):
    dNelect = nelect - val
    wkFunc = -fermi_ene - fermi_shift
    pot_she = wkFunc - u_ref
    # Bugged VASPsol? surrently needed
    energy += fermi_shift * (dNelect)
    ene_corr = energy + dNelect * wkFunc
    return pot_she, ene_corr

def make_surfChrg_sp(nelect):
    homedir = os.getcwd()
    os.system(f'mkdir n_{nelect:.2f}')
    os.chdir(f'n_{nelect:.2f}')
    os.system('cp ../../KPOINTS ../POSCAR ../POTCAR .')
    os.system('cp ../../INCAR-sc INCAR')
    with open('INCAR', 'a') as f:
        f.write(f'NELECT={nelect}')
    os.chdir(homedir)

def make_surfChrg_batch(pp_path, list_deltaCharge, shift=0):
    pos2pot(pp_path)
    nelect_neu = get_neu_nelect(shift=shift)
    for d in list_deltaCharge:
        nelect_tmp = nelect_neu + d
        make_surfChrg_sp(nelect_tmp)

def do_surfChrg_sp(nelect, vasp_cmd, u_ref=4.44, shift=0):
    homedir = os.getcwd()
    os.system(f'mkdir n_{nelect:.2f}')
    os.chdir(f'n_{nelect:.2f}')
    os.system('cp ../../KPOINTS ../POSCAR ../POTCAR .')
    os.system('cp ../../INCAR-sc INCAR')
    with open('INCAR', 'a') as f:
        f.write(f'NELECT={nelect}')
    os.system(vasp_cmd)
    val, nelect, ene, efermi, shftfermi = extractVASPsol(shift=shift)
    USHE, G = pb_calc(val, nelect, ene, efermi, shftfermi, u_ref=u_ref)
    os.system(
        'rm WAVECAR CHG CHGCAR vasprun.xml POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')
    os.chdir(homedir)
    return USHE, G


def do_surfChrg_batch(pp_path, list_deltaCharge, vasp_cmd, u_ref=4.44, shift=0):
    pos2pot(pp_path)
    nelect_neu = get_neu_nelect(shift=shift)
    for d in list_deltaCharge:
        nelect_tmp = nelect_neu + d
        ushe, g = do_surfChrg_sp(nelect_tmp, vasp_cmd, u_ref=u_ref, shift=shift)
        with open('sc.dat', 'a') as f:
            f.write(f'{d}\t{ushe}\t{g}\n')

def do_surfChrg_2step_sp(nelect, vasp_cmd, d, u_ref=4.44):
    homedir = os.getcwd()
    nelect_old = nelect - d
    os.system(f'mkdir n_{nelect:.2f}')
    os.chdir(f'n_{nelect:.2f}')
    os.system('cp ../POSCAR ../POTCAR .')
    if d!=0:
        os.system(f'cp ../n_{nelect_old:.2f}/WAVECAR ../n_{nelect_old:.2f}/CHGCAR .')
    for i in range(1,3):
        os.system('cp ../../INCAR-%i INCAR'%i)
        os.system('cp ../../KPOINTS-%i KPOINTS'%i)
        with open('INCAR', 'a') as f:
            f.write(f'NELECT={nelect}')
        os.system(vasp_cmd)
    val, nelect, ene, efermi, shftfermi = extractVASPsol()
    USHE, G = pb_calc(val, nelect, ene, efermi, shftfermi, u_ref=u_ref)
    os.system(
        'rm CHG vasprun.xml POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')
    os.chdir(homedir)
    return USHE, G

def do_surfChrg_batch_wGuess(pp_path, list_deltaCharge, vasp_cmd, u_ref=4.44):
    pos2pot(pp_path)
    nelect_neu = get_neu_nelect()
    for d in list_deltaCharge:
        nelect_tmp = nelect_neu + d
        ushe, g = do_surfChrg_2step_sp(nelect_tmp, vasp_cmd, d, u_ref=u_ref)
        with open('sc.dat', 'a') as f:
            f.write(f'{d}\t{ushe}\t{g}\n')
    os.system('rm */CHGCAR */WAVECAR')
            
def get_parabola(dir_sc='.'):
    import matplotlib.pyplot as plt
    from sklearn.metrics import r2_score
    data = np.loadtxt(f'{dir_sc}/sc.dat')
    x = data[:, 1]
    y = data[:, 2]
    a, b, c = np.polyfit(x, y, 2)
    r_sqr = r2_score(y, a*x**2 + b*x + c)
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = a*x_fit**2 + b*x_fit + c
    plt.scatter(x, y, label='DFT')
    plt.plot(x_fit, y_fit, label='fit')
    plt.xlabel('U$_\mathrm{SHE}$ (V)', fontsize='x-large')
    plt.ylabel('G (eV)', fontsize='x-large')
    plt.legend()
    textstr = '%.6f, %.6f, %.6f; R^2 = %.6f' % (a, b, c, r_sqr)
    plt.title(textstr)
    fn = os.getcwd().split('/')[-1]
    plt.savefig(f'fit-{fn}.png', bbox_inches="tight")
    plt.close()
    with open(f'{dir_sc}/parabola.dat', 'w') as f:
        f.write(f'{a}\t{b}\t{c}')
