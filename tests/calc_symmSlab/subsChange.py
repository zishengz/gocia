from ase.io import read, write
import sys

mainName = sys.argv[1]
subName = sys.argv[2]
cutZ = eval(sys.argv[3])

sMain = read(mainName)
sSubs = read(subName)

del sMain[[a.index for a in sMain if a.position[2] < cutZ]]
del sSubs[[a.index for a in sSubs if a.position[2] > cutZ]]

sMain.extend(sSubs)

write('subs-'+mainName, sMain)

