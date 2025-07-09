from gocia_new.interface import Interface
from gocia_new.geom import build 
from ase.io import read,write
from ase import Atoms

subAtoms = read('sub-POSCAR')
atoms = read('group.db',':')
atoms = atoms[0]

tmpint = Interface(atoms,subAtoms,fragList=[[193, 195], [194, 196]],clusterAtoms=['Rh','Rh'])
clusterList = tmpint.get_cluster_fragList()
clusterList = [item for row in clusterList for item in row]

clusterAtom = Atoms()
for i in clusterList:
    clusterAtom.extend(atoms[i])

clusterAtom.set_positions(clusterAtom.get_positions() - clusterAtom.get_positions()[0])

newint = build.boxSample_mol(tmpint,[clusterAtom],xyzLims=tmpint.get_sampling_box(),doShuffle=False)
finalcoords = tmpint.get_allAtoms().get_positions()
finalcoords[clusterList,:] =  newint.get_allAtoms().get_positions()[-len(clusterList):]

final = tmpint.get_allAtoms()
final.set_positions(finalcoords)

tmpint.set_allAtoms(final)

