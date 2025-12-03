# Copyright (c) 2023, Zisheng Zhang (SUNCAT Center)
# All rights reserved.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

"""
Description:
    This script generates a solvent layer around a central molecular interface.
    It randomly places solvent molecules and a specified cation within a defined
    simulation box, checking for steric clashes to ensure a physically realistic
    arrangement. The script is designed to generate multiple configurations and
    save them as a trajectory file.

Usage:
    python solvSamp.py <input_xyz_file>

    - <input_xyz_file>: An XYZ file containing the initial interface structure.
"""

import sys
import random
import numpy as np
from scipy.spatial.transform import Rotation as R

from ase import Atoms
from ase.build import molecule
from ase.data import vdw_radii
from ase.io import read, write

from gocia.interface import Interface
from gocia import geom


def get_vdwRadii(atoms):
    """
    Returns a list of van der Waals radii for each atom in an ASE Atoms object.
    Includes a special case for Manganese (25) and Iron (26).
    """
    # Provide a default radius of 0.75 for Mn and Fe.
    return [vdw_radii[i] if i not in [25, 26] else 0.75 for i in atoms.numbers]


def has_intermolecular_clash(mol_a, mol_b, tolerance=0, allowed_bonds=None):
    """
    Checks for steric clashes between two ASE Atoms objects (mol_a and mol_b).

    A clash is defined as the distance between two atoms being smaller than the
    sum of their van der Waals radii, minus a tolerance. A "very bad contact"
    (distance < 0.5 * sum of radii) will always be treated as a clash,
    overriding the `allowed_bonds` list.

    Args:
        mol_a (Atoms): The first molecule.
        mol_b (Atoms): The second molecule, which should contain the periodic cell info.
        tolerance (float): A tolerance value to subtract from the VDW radii sum.
        allowed_bonds (list): A list of atom symbol pairs (e.g., [['Cs', 'O']])
                              that are exempt from normal clash detection.

    Returns:
        bool: True if a clash is detected, False otherwise.
    """
    if not len(mol_a) or not len(mol_b):
        return False
    if allowed_bonds is None:
        allowed_bonds = []

    # Use a set for efficient lookup of allowed bonds.
    allowed_set = {tuple(sorted(p)) for p in allowed_bonds}

    radii_a = get_vdwRadii(mol_a)
    radii_b = get_vdwRadii(mol_b)
    symbols_a = mol_a.get_chemical_symbols()
    symbols_b = mol_b.get_chemical_symbols()

    # To correctly handle Minimum Image Convention (MIC) for intermolecular
    # distances, combine the molecules into a single Atoms object. mol_b,
    # which comes from the main interface, contains the cell definition.
    combined = mol_b + mol_a
    all_dists = combined.get_all_distances(mic=True)

    # Extract the sub-matrix corresponding to distances between mol_a and mol_b.
    dists = all_dists[len(mol_b):, :len(mol_b)]

    # Check each pair of atoms (one from each molecule).
    for i, r_a in enumerate(radii_a):
        for j, r_b in enumerate(radii_b):
            sum_radii = r_a + r_b

            # 1. Check for "very bad contacts" which are always clashes.
            if dists[i, j] < 0.5 * sum_radii:
                return True  # This clash overrides the allowed_bonds list.

            # 2. Check for normal clashes.
            if dists[i, j] < (sum_radii - tolerance):
                # If a potential clash is found, check if it's an allowed bond.
                pair = tuple(sorted((symbols_a[i], symbols_b[j])))
                if pair in allowed_set:
                    continue  # It's an allowed bond, so ignore it.

                return True  # It's a real clash.
    return False


def gen_solv_layer(interfc, solvList, xyzLims, toler=0, allowed_bonds=None):
    """
    Generates a solvent layer by randomly placing molecules from solvList
    around the interface.

    Returns:
        Interface: The solvated interface object, or None if a molecule
                   cannot be placed within the allowed number of attempts.
    """
    tmpInterfc = interfc.copy()
    random.shuffle(solvList)

    n_solvents = len(solvList)
    n_added = 0
    total_trials = 0

    for i in range(n_solvents):
        solv_mol = solvList[i].copy()

        placed = False
        n_trials_for_this_mol = 0
        max_trials_per_mol = 1500  # Avoids getting stuck on a single molecule.

        while not placed:
            if n_trials_for_this_mol >= max_trials_per_mol:
                print(f"\nFailed to place molecule {n_added + 1} after "
                      f"{max_trials_per_mol} attempts. Aborting.")
                return None

            total_trials += 1
            n_trials_for_this_mol += 1

            # 1. Generate a trial position and orientation for the molecule.
            centreCoord = geom.rand_point_box(xyzLims)
            rot_axis = geom.rand_direction()
            rot_angle = random.random() * np.pi
            rot = R.from_rotvec(-rot_angle * rot_axis)

            trial_mol = solv_mol.copy()
            trial_mol.set_positions(rot.apply(trial_mol.get_positions()))
            trial_mol.translate(centreCoord)

            # 2. Check for clashes with all atoms already in the interface.
            existing_atoms = tmpInterfc.get_allAtoms()
            if not has_intermolecular_clash(
                trial_mol, existing_atoms, toler, allowed_bonds=allowed_bonds
            ):
                # 3. If no clash, add the molecule to the interface.
                tmpInterfc.merge_adsorbate(trial_mol)
                n_added += 1
                placed = True
                print(f'{n_added}({total_trials},{n_trials_for_this_mol})', end='  ', flush=True)

    print(f'\nAcceptance rate: {n_solvents / total_trials * 100:.2f}%')
    return tmpInterfc


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python solvSamp.py <input_slab_file> <number of samples>")
        sys.exit(1)

    inpname = sys.argv[1]
    cation = 'Na'

    # Initialize the surface interface.
    surf = Interface(inpname, inpname)
    surf.print()

    # Define atom pairs that can be closer than their VDW radii sum.
    # This is useful for modeling coordination bonds (e.g., cation-oxygen).
    # allowed_bonds = [[cation, 'O'], [cation, 'H']]

    traj = []
    num_configs_to_generate = int(sys.argv[2])

    while len(traj) < num_configs_to_generate:
        print(f'\nGenerating configuration [{len(traj)+1}/{num_configs_to_generate}]', end='\t')

        # Get the lateral dimensions from the simulation cell.
        cell = surf.get_allAtoms().get_cell()

        simulation_box = np.array([
            [0, cell[0][0]],
            [0, cell[1][1]],
            [7, 12]]) # you need to specify the Z range


        # Place the solvent molecules to solvate
        slab_solvated = gen_solv_layer(
            surf,
            [molecule('C3H8')],
            simulation_box,
            toler=0.75,
        )

        if slab_solvated is None:
            print("Restarting configuration generation...")
            continue  # If solvation failed, try again.

        # Add the completed structure to the trajectory.
        slab_solvated.sort()
        traj.append(slab_solvated.get_allAtoms())

    # Write the complete trajectory to a database file.
    if traj:
        print(f"\nSuccessfully generated {len(traj)} configurations.")
        write('ini.db', traj)
        print("Trajectory saved to ini.db")
    else:
        print("\nFailed to generate any valid configurations.")
