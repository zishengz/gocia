
import numpy as np
from ase.io import read, write
from ase.formula import Formula
from ase.units import kB

k_B = 8.61733E-05 #eV/K

class GCE:
    # every Atoms object in here MUST have ase calculator and results

    def __init__(self, traj):
        self.traj = traj

        elems = []
        for a in self.traj:
            elems = set(list(elems) + a.get_chemical_symbols())
        self.elems = list(set(elems))

        elem_counts = {}
        for i in range(len(self.elems)):
            list_tmp = []
            for a in self.traj:
                list_tmp.append(a.get_chemical_symbols().count(self.elems[i]))
            elem_counts[self.elems[i]] = list_tmp
        self.elem_counts = elem_counts

        self.epots = [a.get_potential_energy() for a in self.traj]

    def get_gcfe(self, atoms_sub, dict_mu, indices=None):
        if indices is None:
            indices = range(len(self.epots))
        if type(atoms_sub) is str:
            atoms_sub = read(atoms_sub)
        elems_gc = list(dict_mu)
        list_gcfe = []
        for i in indices:
            tmp = self.epots[i] - atoms_sub.get_potential_energy()
            for e in elems_gc:
                tmp -= (self.elem_counts[e][i]-atoms_sub.get_chemical_symbols().count(e))*dict_mu[e]
            list_gcfe.append(tmp)
        return list_gcfe


    def get_stoi_grps(self):
        grp_stoi = []
        grp_ind = []
        for i in range(len(self.traj)):
            a = self.traj[i]
            stoi = a.get_chemical_formula()
            if stoi not in grp_stoi:
                grp_stoi.append(stoi)
                grp_ind.append([i])
            else:
                grp_ind[grp_stoi.index(stoi)].append(i)

        for i in range(len(grp_stoi)):
            l_ind = grp_ind[i]
            l_ene = [self.epots[m] for m in l_ind]
            _, l_ind = zip(*sorted(zip(l_ene, l_ind)))
            grp_ind[i] = l_ind

        return grp_stoi, grp_ind
    
    def get_GM_index(self, grp=None):
        if grp is None:
            print('MUST HAVE A GROUP (list of structue id\'s) AS INPUT!')
            return None
        else:
            tmp = [self.epots[m] for m in grp]
            return self.epots.index(min(tmp))
        
    def get_pop(self, energy, grp, T=298.15):
        # a single energy value for each structure
        # TODO: X-dependent energy
        if len(energy) != len(grp):
            print('The lengths of energy and indices lists must be consistent!')
        if type(energy) is list:
            energy_zeroed = np.array(energy) - min(energy)
        elif type(energy) is np.ndarray:
            energy_zeroed = energy - energy.min()
        else:
            print('BAD ENERGY! it should be either a list or array')
        
        pop = np.exp(-energy_zeroed / k_B / T)
        pop /= pop.sum(axis=0).T
        return pop
    
    def get_pop_1d(self, energy, grp, T=298.15):
        # a single energy value for each structure
        # TODO: X-dependent energy
        if len(energy) != len(grp):
            print('The lengths of energy and indices lists must be consistent!')
        if type(energy) is np.ndarray:
            energy_zeroed = np.array(energy) - energy.min(axis=0)
        else:
            print('BAD ENERGY! it should be either a 2d array')
        
        pop = np.exp(-energy_zeroed / k_B / T)
        pop /= pop.sum(axis=0).T
        return pop
    
    




        

#    def get LELM and sort

            



# # def get_GM(traj):
    
    

# # def get_LELM():
    

# # def get_pop():
    

# # def get_ens_avg():


