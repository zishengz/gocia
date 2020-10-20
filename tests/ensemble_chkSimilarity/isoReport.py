
import os, sys
from ase.db import connect
from ase.io import read
from gocia.interface import Interface
import gocia.utils.report as rp
from gocia.utils.visualize import histogram

substrate = read(sys.argv[1])
srtDB = connect(sys.argv[2])

baseName = ''
eneList = []
for r in srtDB.select([('eV2GM','<', 5.0)]):
    a = r.toatoms()
    eneList.append(r.eV)
    surf = Interface(
        allAtoms=a,
        subAtoms=substrate,
    )
    baseName = surf.tags
    surf.tags += '-%i'%r.id
    surf.draw(
        'BS',
        title='#%i, E = %.3f eV, mag = %.2f'%\
            (r.id, r.eV2GM, r.mag),
    )

# myImgList = [f for f in os.listdir() if baseName in f and 'png' in f]
# rp.gen_pdf_test(
#     projName='Global Optimization of Restructured %s'%baseName,
#     description=\
#         '\\begin{itemize}\n'+\
#         '\\item From direct BLDA sampling with random rattling of 0.1 Angstrom and Hookean pre-optimization.\n'+\
#         '\\item The same-element probability penalty is set to 50\%, and desorbed N$_2$ structures are filtered out.\n'+\
#         '\\item Color code: Ga - light brown; N - blue.\n'+\
#         '\\item For clarity, the constrained atoms are drawn in white, and the buffer atoms are drawn in shallow color.\n'+\
#         '\\end{itemize}',
#     imgList=myImgList,
# )
# os.system('rm *%s*.png'%baseName)

from natsort import natsorted
myImgList = natsorted([f for f in os.listdir() if baseName in f and 'png' in f and 'bs' in f])
# histogram(eneList, baseName)
myBigFig = natsorted([f for f in os.listdir() if baseName in f and 'png' in f and 'hist' in f])

# # rp.gen_pdf_test2(
# #     projName='Global Optimization of Restructured %s'%baseName,
# #     description=\
# #         '\\begin{enumerate}\n'+\
# #         '\\item From direct BLDA sampling with random rattling of 0.1 Angstrom and Hookean pre-optimization.\n'+\
# #         '\\item The same-element probability penalty is set to 50\%, and desorbed N$_2$ structures are filtered out.\n'+\
# #         '\\item Color code: Ga - light brown; N - blue.\n'+\
# #         '\\item For clarity, the constrained atoms are drawn in white, and the buffer atoms are drawn in shallow color.\n'+\
# #         '\\end{enumerate}',
# #     imgList=myImgList,
# #     bigFig=myBigFig
# # )

rp.gen_pdf_GlobOpt1(
    projName='Global Optimization of Restructured %s'%baseName,
    description='comment.tex',
    imgList=myImgList,
    bigFig=myBigFig
)

os.system('rm *%s*.png'%baseName)