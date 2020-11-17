from gocia.interface import Interface
from gocia.data import covalRadii
from ase.io import read, write
from ase.atoms import Atoms
import numpy as np

def get_surf_byZ(atoms, zRange=[0, 999]):
    myPos = atoms.get_positions()
    return [
        i for i in range(len(atoms))\
        if min(zRange) < myPos[i][2] < max(zRange)
    ]

def get_surfAtom_byXY(atoms, xypos, scale=2):
    radiiAll = scale * np.array([
        covalRadii[i] for i in atoms.get_atomic_numbers()])
    z_test = max(atoms.get_positions()[:, 2])+5
    surfIndex = []
    while len(surfIndex) == 0:
        testAtom = Atoms(['H'], [[xypos[0], xypos[1], z_test]])
        tmpAtoms = atoms.copy()
        tmpAtoms.extend(testAtom)
        dist_test = list(tmpAtoms.get_distances(-1, list(range(len(atoms)))) - radiiAll)
        for i in dist_test:
            if i < 0: surfIndex.append(dist_test.index(i))
        z_test -= 0.05
    return surfIndex

def get_surf_grid(atoms, mesh=10):
    cell = atoms.get_cell()
    cellScale = np.linspace(0, 1, mesh, endpoint=False)
    surfList = []
    for i,j in zip(cellScale, cellScale):
        surfList += get_surfAtom_byXY(atoms,\
            [cell[0][0]*i+cell[1][0]*j, cell[0][2]*i+cell[1][1]*j])
    # for i in cellScale:
    #     for j in cellScale:
    #         surfList += get_surfAtom_byXY(atoms,\
    #             [cell[0][0]*i+cell[1][0]*j, cell[0][2]*i+cell[1][1]*j])
    return sorted(list(set(surfList)))

def get_extended_atoms(atoms):
    tmpAtoms = atoms.copy()
    # tmpAtoms.set_positions(tmpAtoms.get_positions()\
    #     -atoms.get_cell()[0]/2-atoms.get_cell()[1]/2)
    tmpAtoms = atoms*[3,3,1]
    tmpAtoms.set_positions(tmpAtoms.get_positions()\
        -atoms.get_cell()[0]/1-atoms.get_cell()[1]/1)
    return tmpAtoms

# TODO DUPLICATE CHECK
# TODO ADD ATOP
# TODO ADD HOLLOW
# TODO ADD BRIDGE

surf = Interface(
    allAtoms='calcSlab.vasp',
    subAtoms='sub-calcSlab.vasp'
)

surf.print()

surf = read('out-Cu111.vasp')
#surf = read('WBsurf2.vasp')

#print(get_surf_byZ(surf, [9.5, 10.5]))
#print(get_surf_grid(surf, 20))

surfList = get_surf_byZ(surf, [16, 18])
#surfList = get_surf_grid(surf, 20)
print(surfList)
del surf[[i for i in range(len(surf)) if i not in surfList]]
surf_ext = get_extended_atoms(surf)
#surf_ext = surf

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
from scipy.spatial import Delaunay

pos_ext = surf_ext.get_positions()
tri = Delaunay(pos_ext[:, :2])
tri_nodes = tri.simplices
pos_nodes = pos_ext[tri_nodes]
hollow = np.array([t.sum(axis=0)/3 for t in pos_nodes])
print(hollow)
bridge = []
for i in pos_nodes:
    bridge.append((i[0]+i[1])/2)
    bridge.append((i[0]+i[2])/2)
    bridge.append((i[1]+i[2])/2)
bridge = np.array(bridge)

atop = surf.get_positions()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.triplot(pos_ext[:,0], pos_ext[:,1], triangles=tri.simplices, color='grey',)
ax.plot(hollow[:,0],  hollow[:,1],  'ok', label='Hollow')
ax.plot(bridge[:,0],  bridge[:,1],  'or', label='Bridge')
ax.plot(pos_ext[:,0], pos_ext[:,1], 'ob', label='Atop')

plt.axis('scaled')
cell_param = surf.get_cell()
ax.set_xlim([cell_param[1][0],cell_param[0][0]])
ax.set_ylim([cell_param[0][1],cell_param[1][1]])

plt.legend(loc='lower left')
plt.savefig('ads.png',  bbox_inches = "tight", transparent=True)

ads_coord = np.array([[0, 0, 0], [0, 0, 1.17]])
#surf.extend(Atoms('CO', ads_coord + hollow[0] + np.array([0,0,1.5])))
#surf.extend(Atoms('CO', ads_coord + bridge[0] + np.array([0,0,1.5])))
#surf.extend(Atoms('CO', ads_coord + atop[0] + np.array([0,0,1.5])))
surf.wrap()
write('surf.png', surf)
write('surf.vasp', surf)

# TODO find ads direction from normal vector


# TODO wrap back to original cell from ext cell
