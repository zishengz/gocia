from __future__ import division

import numpy as np

from gocia.data import covalRadii
from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError


class LennardJones(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']
    default_parameters = {'epsilon': 0.1,
                          'sigma': 1.0,
                          'tolerAngs': 0.2,
                          'tolerMult': 0.2,
                          'cutoff': None}
    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        sigma = self.parameters.sigma
        epsilon = self.parameters.epsilon
        tolerAngs = self.parameters.tolerAngs
        tolerMult = self.parameters.tolerMult        
        cutoff = self.parameters.cutoff
        if cutoff is None:
            cutoff = 0.9

        if 'numbers' in system_changes:
            self.nl = NeighborList([cutoff / 2] * natoms,\
                      self_interaction=False)

        self.nl.update(self.atoms)

        positions = self.atoms.positions
        cell = self.atoms.cell

        e0 = 4 * epsilon * ((sigma / cutoff)**12 - (sigma / cutoff)**6)

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

        def getLJparam(elem1, elem2):
            sigma_tmp = covalRadii[elem1] + covalRadii[elem1]
            return sigma_tmp


        for a1 in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            r2 = (d**2).sum(1)
            covalBl = covalRadii[self.atoms.get_atomic_numbers()[a1]]+\
                covalRadii[self.atoms.get_atomic_numbers()[neighbors]]
            c6 = (sigma**2 / r2)**3
            c6[covalBl - np.sqrt(r2) < tolerAngs] = 0.0
            c6[np.sqrt(r2) / covalBl > 1 - tolerMult] = 0.0
#            c6[r2 > radius**2] = 0.0
            energy -= e0 * (c6 != 0.0).sum()
            c12 = c6**2
            energy += 4 * epsilon * (c12 - c6).sum()
            f = (24 * epsilon * (2 * c12 - c6) / r2)[:, np.newaxis] * d
            forces[a1] -= f.sum(axis=0)
            for a2, f2 in zip(neighbors, f):
                forces[a2] += f2
            stress += np.dot(f.T, d)

        if 'stress' in properties:
            if self.atoms.number_of_lattice_vectors == 3:
                stress += stress.T.copy()
                stress *= -0.5 / self.atoms.get_volume()
                self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]
            else:
                raise PropertyNotImplementedError

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces


