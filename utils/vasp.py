import os
import numpy as np
from ase.io import read, write
from gocia.geom import get_fragments, del_freeMol


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


def is_vaspSuccess(jobdir='.'):
    return 'E0' in open(f'{jobdir}/OSZICAR').readlines()[-1]


def do_multiStep_opt(step=3, vasp_cmd='', chkMol=True):
    for i in range(1, step+1):
        os.system(f'cp ../INCAR-{i} INCAR')
        os.system(vasp_cmd)
        if not is_vaspSuccess():
            os.system('touch FAIL')
            exit()
        os.system(f'cp CONTCAR out-{i}.vasp')
        if not chkMol:
            os.system('cp CONTCAR POSCAR')
        else:
            write('POSCAR', del_freeMol(read('CONTCAR')))
            # Make sure the final structure has no free molecule
            if i == step:
                if len(get_fragments(read('CONTCAR'))) > 1:
                    os.system('touch BADSTRUCTURE')
                    exit()
    os.system('rm WAVECAR CHG CHGCAR vasprun.xml POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')


# Below are for surface charging calculations

def get_neu_nelect(poscar='POSCAR', potcar='POTCAR'):
    elem_list = [l.split()[-2]
                 for l in open(potcar).readlines() if 'TITEL' in l]
    elem_list = [e.split('_')[0] if '_' in e else e for e in elem_list]
    zval_list = [eval(l.split()[-4])
                 for l in open(potcar).readlines() if 'ZVAL' in l]
    atom_list = read(poscar).get_chemical_symbols()
    return sum([zval_list[elem_list.index(a)] for a in atom_list])


def extractVASPsol(dirName='.'):
    tmpHome = os.getcwd()
    os.chdir(dirName)
    val = get_neu_nelect()
    fermi_shift = eval([l for l in open('out').readlines()
                        if 'FERMI_SHIFT' in l][-1].split()[2])
    fermi_ene = eval([l for l in open('OUTCAR').readlines()
                      if 'E-fermi' in l][-1].split()[2])
    nelect = eval([l for l in open('INCAR').readlines()
                   if 'NELECT' in l][-1].split('=')[-1])
    energy = eval(open('OSZICAR').readlines()[-1].split()[4])
    os.chdir(tmpHome)
    return val, nelect, energy, fermi_ene, fermi_shift


def pb_calc(val, nelect, energy, fermi_ene, fermi_shift):
    dNelect = nelect - val
    wkFunc = -fermi_ene - fermi_shift
    pot_she = wkFunc - 4.44
    # Bugged VASPsol? surrently needed
    energy += fermi_shift * (dNelect)
    ene_corr = energy + dNelect * wkFunc
    return pot_she, ene_corr


def do_surfChrg_sp(nelect, vasp_cmd):
    homedir = os.getcwd()
    os.system(f'mkdir n_{nelect:.2f}')
    os.chdir(f'n_{nelect:.2f}')
    os.system('cp ../../KPOINTS ../POSCAR ../POTCAR .')
    os.system('cp ../../INCAR-sc INCAR')
    with open('INCAR', 'a') as f:
        f.write(f'NELECT={nelect}')
    os.system(vasp_cmd)
    val, nelect, ene, efermi, shftfermi = extractVASPsol()
    USHE, G = pb_calc(val, nelect, ene, efermi, shftfermi)
    os.system(
        'rm WAVECAR CHG CHGCAR vasprun.xml POTCAR PCDAT XDATCAR DOSCAR EIGENVAL IBZKPT')
    os.chdir(homedir)
    return USHE, G


def do_surfChrg_batch(list_deltaCharge, vasp_cmd):
    pos2pot(input.pp_path)
    nelect_neu = get_neu_nelect()
    for d in list_deltaCharge:
        nelect_tmp = nelect_neu + d
        ushe, g = do_surfChrg_sp(nelect_tmp, vasp_cmd)
        with open('sc.dat', 'a') as f:
            f.write(f'{d}\t{ushe}\t{g}\n')


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
