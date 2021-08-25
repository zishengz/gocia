
from ase.constraints import FixAtoms
from ase.io import read, write
import sys, os

inp = sys.argv[1]
prj = os.getcwd().split('/')[-1]

traj = read(inp, index=':')
for i in range(len(traj)):
#	curName = '%s-%s.vasp'%(prj, str(i).zfill(3))
	curName = 's%s.vasp'%(str(i).zfill(3))
	write(curName, traj[i], vasp5=True)
	print(' > %s written!'%curName)
