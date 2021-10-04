
import numpy as np
from ase.db import connect
from gocia.interface import Interface
from gocia.ga import (
    get_eneFactor,
    get_matedFactor,
    get_matedFactor2,
    get_matedFactor3
)
from gocia.ga.crossover import crossover_snsSurf_2d_GC
from gocia.ensemble.comparator import srtDist_similar_zz
from ase.io import read, write


class PopulationGrandCanonical:
    def __init__(
        self,
        substrate=None,
        gadb=None,
        popSize=20,
        zLim=None,
        subsPot=0,
        chemPotDict=None,
        compParam=None,
        matingMethod=None,
        convergeCrit=None,
        simParam1=5e-4,
        simParam2=0.25
    ):

        if gadb is not None:
            self.gadb = connect(gadb)

        if type(substrate) is str:
            self.substrate = read(substrate)
        else:
            self.substrate = substrate

        if zLim is not None:
            self.zLim = zLim
        if zLim is None:
            self.zLim = Interface(self.substrate, self.substrate).zLim

        self.popSize = popSize

        if convergeCrit is None:
            self.convergeCrit = 5 * self.popSize
        else:
            self.convergeCrit = convergeCrit

        self.subsPot = subsPot

        if chemPotDict is not None:
            self.chemPotDict = chemPotDict

        self.iniSize = len(self)

        self.simParam1 = simParam1
        self.simParam2 = simParam2

    def __len__(self):
        return len(self.gadb)

    def is_converged(self):
        gmid = self.get_GMrow().id
        if gmid < self.iniSize:
            return False
        else:
            return len(self) - gmid > self.convergeCrit

    def is_converged2(self):
        gmid = eval(open('gmid').read())
        if gmid < self.iniSize:
            return False
        else:
            return len(self) - gmid > self.convergeCrit

    def calc_grandPot(self, atoms, dftene):
        #        myRow = self.gadb.get(id=myID)
        #        myPot = myRow['eV']
        myPot = dftene - self.subsPot
        adsSymbol = Interface(atoms, self.substrate).\
            get_adsAtoms().get_chemical_symbols()
        for s in adsSymbol:
            if s in self.chemPotDict:
                myPot -= self.chemPotDict[s]
        return myPot

    def initializeDB(self, sc=False):
        if sc:
            for i in range(len(self.gadb)):
                r = self.gadb.get(id=i+1)
                self.gadb.update(i+1, mated=0, alive=1,
                                grandPot=self.calc_grandPot(r.toatoms(), r.sc_eV/2), label='0 0 init')
        else:
            for i in range(len(self.gadb)):
                r = self.gadb.get(id=i+1)
                self.gadb.update(i+1, mated=0, alive=1,
                                grandPot=self.calc_grandPot(r.toatoms(), r.eV), label='0 0 init')
        self.natural_selection()
        gm = self.get_GMrow()
        with open('gmid', 'w') as f:
            f.write(str(gm.id))

    def get_ID(self, condString):
        tmp = []
        for r in self.gadb.select(condString):
            tmp.append(r.id)
        return tmp

    def get_valueOf(self, valueString, idList):
        tmp = []
        for i in idList:
            tmp.append(
                self.gadb.get(id=i)[valueString]
            )
        return tmp

    def get_GMrow(self):
        eneList = self.get_valueOf('grandPot', self.get_ID('done=1'))
        self.Emin = min(eneList)
        row_min = None
        for r in self.gadb.select(grandPot=self.Emin):
            row_min = r
        return row_min

    def get_fitness(self, idList):
        tmpList = self.get_valueOf('grandPot', idList)
        f_ene = get_eneFactor(tmpList)
        tmpList = self.get_valueOf('mated', idList)
        f_mat = get_matedFactor3(tmpList)
        return f_ene * f_mat

    def choose_parents(self):
        aliveList = self.get_ID('alive=1')
        fitList = self.get_fitness(aliveList)
        weights = fitList / sum(fitList)
        parents = np.random.choice(aliveList, size=2, replace=False, p=weights)
        parents = [int(parents[0]), int(parents[1])]
        return parents

    def natural_selection(self):
        aliveList = self.get_ID('alive=1')
        if len(aliveList) > self.popSize:
            fitnessList = self.get_fitness(aliveList)
            aliveList = [x for _, x in sorted(zip(fitnessList, aliveList))]
            deadList = aliveList[:-self.popSize]
            for d in deadList:
                print('REST IN PEACE, %i!' % d)
                self.gadb.update(d, alive=0)

    def is_uniqueInPop(self, atoms, grandPot, eneCut=0.05):
        '''
        Check similarity against the current population
        '''
        aliveList = self.get_ID('alive=1')
        grandPotList = self.get_valueOf('grandPot', aliveList)
        isUnique = True
        for ii in range(len(aliveList)):
            a = self.gadb.get(id=aliveList[ii]).toatoms()
            if a.get_chemical_formula() == atoms.get_chemical_formula():
                if srtDist_similar_zz(atoms, a, self.simParam1, self.simParam2):
                    if eneCut > 0:
                        if -eneCut < grandPotList[ii] - grandPot < eneCut:
                            isUnique = False
                            break
                    else:
                        isUnique = False
                        break
        return isUnique

    def is_uniqueInAll(self, atoms, grandPot, eneCut=0.05):
        '''
        Check similarity against the all sampled structures
        '''
        aliveList = self.get_ID('done=1')
        grandPotList = self.get_valueOf('grandPot', aliveList)
        isUnique = True
        for ii in range(len(aliveList)):
            a = self.gadb.get(id=aliveList[ii]).toatoms()
            if a.get_chemical_formula() == atoms.get_chemical_formula():
                if srtDist_similar_zz(atoms, a, self.simParam1, self.simParam2):
                    if eneCut > 0:
                        if -eneCut < grandPotList[ii] - grandPot < eneCut:
                            isUnique = False
                            break
                    else:
                        isUnique = False
                        break
        return isUnique

    def is_uniqueInPop_geom(self, atoms):
        '''
        Check similarity against the all sampled structures
        '''
        aliveList = self.get_ID('alive=1')
        isUnique = True
        for ii in range(len(aliveList)):
            a = self.gadb.get(id=aliveList[ii]).toatoms()
            if a.get_chemical_formula() == atoms.get_chemical_formula():
                if srtDist_similar_zz(atoms, a, self.simParam1, self.simParam2):
                    isUnique = False
                    break
        return isUnique

    def is_uniqueInAll_geom(self, atoms):
        '''
        Check similarity against the all sampled structures
        '''
        aliveList = self.get_ID('done=1')
        isUnique = True
        for ii in range(len(aliveList)):
            a = self.gadb.get(id=aliveList[ii]).toatoms()
            if a.get_chemical_formula() == atoms.get_chemical_formula():
                if srtDist_similar_zz(atoms, a, self.simParam1, self.simParam2):
                    isUnique = False
                    break
        return isUnique

    def gen_offspring(self, mutRate=0.3, rattleOn=True, growOn=True, leachOn=True, permuteOn=True, transOn=True, transVec=[[-2, 2], [-2, 2]]):
        kid, parent = None, None
        mater, pater = 0, 0
        while kid is None:
            mater, pater = self.choose_parents()
            a1 = self.gadb.get(id=mater).toatoms()
            a2 = self.gadb.get(id=pater).toatoms()
            surf1 = Interface(a1, self.substrate, zLim=self.zLim)
            surf2 = Interface(a2, self.substrate, zLim=self.zLim)
            kid = crossover_snsSurf_2d_GC(surf1, surf2, tolerance=0.75)
            parent = surf1.copy()
        print('PARENTS: %i and %i' % (mater, pater))
        myMutate = ''
        if srtDist_similar_zz(a1, a2) or not self.is_uniqueInAll_geom(kid.get_allAtoms()):
            print(' |- TOO SIMILAR!')
            mutRate = 1
        if np.random.rand() < mutRate:
            mutType = np.random.choice([0, 1, 2, 3, 4], size=1)[0]
            if mutType == 0 and rattleOn:
                myMutate = 'rattle'
                kid.rattleMut()
            if mutType == 1 and growOn:
                myMutate = 'grow'
                kid.growMut([l for l in self.chemPotDict])
            if mutType == 2 and leachOn:
                myMutate = 'leach'
                kid.leachMut([l for l in self.chemPotDict])
            if mutType == 3 and permuteOn:
                myMutate = 'permute'
                kid.permuteMut()
            if mutType == 4 and transOn:
                myMutate = 'translate'
                kid.transMut(transVec=transVec)
        if len(kid.get_adsList()) <= 1:
            myMutate = 'init'
            print(' |- Bare substrate, BAD!')
            kid = parent.copy()
            kid.rattleMut()
            kid.growMut([l for l in self.chemPotDict])
        open('label', 'w').write('%i %i %s' % (mater, pater, myMutate))
        self.gadb.update(mater, mated=self.gadb.get(id=mater).mated+1)
        self.gadb.update(pater, mated=self.gadb.get(id=pater).mated+1)
        return kid

    def gen_offspring_box(self, mutRate=0.3, xyzLims=[], bondRejList=None, constrainTop=False, rattleOn=True, growOn=True, leachOn=True, permuteOn=True, transOn=True, transVec=[[-2, 2], [-2, 2]]):
        kid, parent = None, None
        mater, pater = 0, 0
        while kid is None:
            mater, pater = self.choose_parents()
            a1 = self.gadb.get(id=mater).toatoms()
            a2 = self.gadb.get(id=pater).toatoms()
            surf1 = Interface(a1, self.substrate, zLim=self.zLim)
            surf2 = Interface(a2, self.substrate, zLim=self.zLim)
            kid = crossover_snsSurf_2d_GC(surf1, surf2, tolerance=0.75)
            parent = surf1.copy()
        print('PARENTS: %i and %i' % (mater, pater))
        myMutate = ''
        if srtDist_similar_zz(a1, a2)\
                or srtDist_similar_zz(a1, kid.get_allAtoms())\
                or srtDist_similar_zz(a2, kid.get_allAtoms()):
            print(' |- TOO SIMILAR!')
            mutRate = 1
        if np.random.rand() < mutRate:
            mutType = np.random.choice([0, 1, 2, 3, 4], size=1)[0]
            if mutType == 0 and rattleOn:
                myMutate = 'rattle'
                kid.rattleMut()
            if mutType == 1 and growOn:
                myMutate = 'grow'
                tmpKid = None
                while tmpKid is None:
                    tmpKid = kid.copy()
                    tmpKid.growMut_box([l for l in self.chemPotDict], xyzLims=xyzLims,
                                bondRejList=bondRejList, constrainTop=constrainTop)
                kid = tmpKid.copy()
            if mutType == 2 and leachOn:
                myMutate = 'leach'
                kid.leachMut([l for l in self.chemPotDict])
            if mutType == 3 and permuteOn:
                myMutate = 'permute'
                kid.permuteMut()
            if mutType == 4 and transOn:
                myMutate = 'translate'
                kid.transMut(transVec=transVec)
        if len(kid.get_adsList()) <= 1:
            myMutate = 'init'
            print(' |- Bare substrate, BAD!')
            kid = parent.copy()
            kid.rattleMut()
            kid.growMut([l for l in self.chemPotDict])
        open('label', 'w').write('%i %i %s' % (mater, pater, myMutate))
        self.gadb.update(mater, mated=self.gadb.get(id=mater).mated+1)
        self.gadb.update(pater, mated=self.gadb.get(id=pater).mated+1)
        return kid

    def add_vaspResult(self, vaspdir='.'):
        import os
        cwdFiles = os.listdir(vaspdir)
        if 'OSZICAR' in cwdFiles\
                and 'BADSTRUCTURE' not in cwdFiles\
                and 'ERROR' not in open('%s/OSZICAR' % vaspdir).read():
            if 'E0' in open(vaspdir+'/OSZICAR', 'r').readlines()[-1]:
                s = read('%s/OUTCAR' % vaspdir, index='-1')
                dirname = os.getcwd().split('/')[-1]
                try:
                    mag = s.get_magnetic_moment()
                except:
                    mag = 0
                ene_eV = s.get_potential_energy()
                grndPot = self.calc_grandPot(s, ene_eV)
                myLabel = open('label', 'r').read()
                print('\n%s IS BORN with G = %.3f eV\t[%s]' % (
                    dirname, grndPot, myLabel))
                if self.is_uniqueInAll(s, grndPot):
                    if grndPot < self.get_GMrow()['grandPot']:
                        print(f' |- {dirname} is the new GM!')
                        with open('gmid', 'w') as f:
                            f.write(str(len(self)))
                    self.gadb.write(
                        s,
                        name=dirname,
                        mag=mag,
                        eV=ene_eV,
                        grandPot=grndPot,
                        mated=0,
                        done=1,
                        alive=1,
                        label=myLabel
                    )
                else:
                    print(f' |- {dirname} is a duplicate!')
                    self.gadb.write(
                        s,
                        name=dirname,
                        mag=mag,
                        eV=ene_eV,
                        grandPot=grndPot,
                        mated=0,
                        done=1,
                        alive=0,
                        label=myLabel
                    )


    def add_vaspResult_SC(self, u_she=0, vaspdir='.'):
        import os
        cwdFiles = os.listdir(vaspdir)
        if 'parabola.dat' in cwdFiles:
            a, b, c = np.loadtxt('parabola.dat')
            ene_sc = a*u_she**2 + b*u_she + c
            # below are non-sc results
            s = read('%s/OUTCAR' % vaspdir, index='-1')
            dirname = os.getcwd().split('/')[-1]
            grndPot = self.calc_grandPot(s, ene_sc/2)
            try:
                mag = s.get_magnetic_moment()
            except:
                mag = 0
            ene_eV = s.get_potential_energy()
            myLabel = open('label', 'r').read()
            print('\n%s IS BORN with G = %.3f eV\t[%s]' % (
                dirname, grndPot, myLabel))
            if self.is_uniqueInAll(s, grndPot):
                if grndPot < self.get_GMrow()['grandPot']:
                    print(f' |- {dirname} is the new GM!')
                    with open('gmid', 'w') as f:
                        f.write(str(len(self)))
                self.gadb.write(
                    s,
                    name=dirname,
                    mag=mag,
                    eV=ene_eV,
                    sc_U=u_she,
                    sc_eV=ene_sc,
                    grandPot=grndPot,
                    a=a,
                    b=b,
                    c=c,
                    mated=0,
                    done=1,
                    alive=1,
                    label=myLabel
                )
            else:
                print(f' |- {dirname} is a duplicate!')
                self.gadb.write(
                    s,
                    name=dirname,
                    mag=mag,
                    eV=ene_eV,
                    sc_U=u_she,
                    sc_eV=ene_sc,
                    grandPot=grndPot,
                    a=a,
                    b=b,
                    c=c,
                    mated=0,
                    done=1,
                    alive=0,
                    label=myLabel
                )

# TODO convergence: TEST
# TODO XTB interface
# TODO CP2K interface
