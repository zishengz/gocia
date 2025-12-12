from ase.io import read
from ase.db import connect
from gocia_new.interface import Interface

subAtoms = read('sub-POSCAR')

db = connect('ini.db')
row = db.get(id=15)
a = row.toatoms()
fragList = row.adsFrags

tmpint = Interface(a,subAtoms,fragList=fragList,clusterAtoms=['Rh','Rh'])

print(tmpint.clusterAtoms)
print()
print(tmpint.get_clusterList())

print(a.get_chemical_symbols()[177])
print(a.get_chemical_symbols()[178])

print(tmpint.get_cluster_fragList(scale=1.2))

print(tmpint.get_fragList())

