from __future__ import division

import numpy as np

from gocia.data import covalRadii
from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError


class Hooke(Calculator):
    implemented_properties = ['energy', 'forces']
    default_parameters = {
                          'k': 1,
                          'tolerAngs': 0.2,
                          'tolerMult': 0.2,
                          'cutoff': 1.5,
                          }
    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        k = self.parameters.k
        tolerAngs = self.parameters.tolerAngs
        tolerMult = self.parameters.tolerMult        
        cutoff = self.parameters.cutoff

        if 'numbers' in system_changes:
            # Want agregation?  cutoff > 2.00
            # Want uniform?     cutoff < 1.50
            # Directional?      cutoff < 1.25
            self.nl = NeighborList([cutoff / 2] * natoms,\
                      self_interaction=False)

        self.nl.update(self.atoms)

        positions = self.atoms.positions
        cell = self.atoms.cell

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

        for a1 in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(a1)
#            print(neighbors)
            cells = np.dot(offsets, cell)
            vecBond = positions[neighbors] + cells - positions[a1]
            lenBond = np.sqrt((vecBond**2).sum(1))
            covalBl = covalRadii[self.atoms.get_atomic_numbers()[a1]]+\
                covalRadii[self.atoms.get_atomic_numbers()[neighbors]]
            compress = covalBl - lenBond
            energy += k / 2 * ((compress)**2).sum()
            # only repulsion
            compress[covalBl - lenBond < tolerAngs] = 0.0
            compress[lenBond / covalBl > 1 -tolerMult] = 0.0
            # # also attraction
            # compress[np.abs(covalBl - lenBond) < tolerAngs] = 0.0
            # compress[np.abs(1 - lenBond / covalBl) < tolerMult] = 0.0
            f = k * (compress)[:, np.newaxis] * vecBond
            forces[a1] -= f.sum(axis=0)
            for a2, f2 in zip(neighbors, f):
                forces[a2] += f2

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces
