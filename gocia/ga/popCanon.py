
import numpy as np
from ase.db import connect
from gocia.interface import Interface
from gocia.ga import (
    get_eneFactor,
    get_matedFactor,
    get_matedFactor2,
    get_matedFactor3
)
from gocia.ga.crossover import crossover_snsSurf_2d
from gocia.ensemble.comparator import srtDist_similar_zz
from ase.io import read, write


class PopulationCanonical:
    def __init__(
        self,
        substrate=None,
        gadb=None,
        popSize=20,
        zLim=None,
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

    def initializeDB(self):
        for i in range(len(self.gadb)):
            self.gadb.update(i+1, mated=0, alive=1)
        self.natural_selection()

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
        eneList = self.get_valueOf('eV', self.get_ID('done=1'))
        self.Emin = min(eneList)
        return self.gadb.get(eV=self.Emin)

    def get_fitness(self, idList):
        tmpList = self.get_valueOf('eV', idList)
        f_ene = get_eneFactor(tmpList)
        tmpList = self.get_valueOf('mated', idList)
        f_mat = get_matedFactor(tmpList)
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

    def is_uniqueInPop(self, atoms, ene):
        '''
        Check similarity against the current population
        '''
        eneCut = 0.05
        aliveList = self.get_ID('alive=1')
        eneList = self.get_valueOf('eV', aliveList)
        isUnique = True
        for ii in range(len(aliveList)):
            a = self.gadb.get(id=aliveList[ii]).toatoms()
            if srtDist_similar_zz(atoms, a, self.simParam1, self.simParam2):
                if -eneCut < eneList[ii] - ene < eneCut:
                    isUnique = False
                    break
        return isUnique

    def gen_offspring(self, mutRate=0.4, rattleOn=True, zEnhance=True, transOn=True, transVec=[[-2, 2], [-2, 2]]):
        kid = None
        mater, pater = 0, 0
        while kid is None:
            mater, pater = self.choose_parents()
            a1 = self.gadb.get(id=mater).toatoms()
            a2 = self.gadb.get(id=pater).toatoms()
            surf1 = Interface(a1, self.substrate, zLim=self.zLim)
            surf2 = Interface(a2, self.substrate, zLim=self.zLim)
            kid = crossover_snsSurf_2d(surf1, surf2, tolerance=0.75)
        print('PARENTS: %i and %i' % (mater, pater))
        if srtDist_similar_zz(a1, a2)\
                or srtDist_similar_zz(a1, kid.get_allAtoms())\
                or srtDist_similar_zz(a2, kid.get_allAtoms()):
            print(' |- TOO SIMILAR!')
            mutRate = 1
        myMutate = ''
        if np.random.rand() < mutRate:
            mutType = np.random.choice([0, 1], size=1)[0]
            if mutType == 0 and rattleOn:
                myMutate = 'rattle'
                kid.rattleMut(zEnhance=zEnhance)
            if mutType == 1 and transOn:
                myMutate = 'translate'
                kid.transMut(transVec=transVec)
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
                mag = s.get_magnetic_moment()
                ene_eV = s.get_potential_energy()
                print('\nA CHILD IS BORN with G = %.3f eV' % (ene_eV))
                if self.is_uniqueInPop(s, ene_eV):
                    if ene_eV < self.get_GMrow()['eV']:
                        print(' |- it is the new GM!')
                    self.gadb.write(
                        s,
                        name=dirname,
                        mag=mag,
                        eV=ene_eV,
                        mated=0,
                        done=1,
                        alive=1
                    )
                else:
                    print(' |- it is a duplicate!')
                    self.gadb.write(
                        s,
                        name=dirname,
                        mag=mag,
                        eV=ene_eV,
                        mated=0,
                        done=1,
                        alive=0
                    )

    def add_lmpResult(self, lmpdir='.'):
        import os
        import gocia.utils.lammps as lmp
        cwdFiles = os.listdir(lmpdir)
        if 'lmp.out' in cwdFiles and 'BADSTRUCTURE' not in cwdFiles:
            if 'Final' in open(lmpdir+'/lmp.out').read():
                dirname = os.getcwd().split('/')[-1]
                mag = 0
                ene_eV = lmp.get_ene(lmpdir+'/lmp.out')
                s = lmp.get_last_frame(lmpdir+'/traj.xyz', lmpdir+'/inp.vasp')
                s.wrap()
                print('\nA CHILD IS BORN with G = %.3f eV' % (ene_eV))
                if self.is_uniqueInPop(s, ene_eV):
                    if ene_eV < self.get_GMrow()['eV']:
                        print(' |- it is the new GM!')
                    self.gadb.write(
                        s,
                        name=dirname,
                        mag=mag,
                        eV=ene_eV,
                        mated=0,
                        done=1,
                        alive=1
                    )
                else:
                    print(' |- it is a duplicate!')
                    self.gadb.write(
                        s,
                        name=dirname,
                        mag=mag,
                        eV=ene_eV,
                        mated=0,
                        done=1,
                        alive=0
                    )


    def add_aseResult(self, atoms, workdir='.'):
        ene_eV = atoms.get_potential_energy()
        print('\nA CHILD IS BORN with G = %.3f eV' % (ene_eV))
        if self.is_uniqueInPop(atoms, ene_eV):
            if ene_eV < self.get_GMrow()['eV']:
                print(f' |- it is the new GM!')
            self.gadb.write(
                atoms,
                name=workdir,
                mag=0,
                eV=ene_eV,
                mated=0,
                done=1,
                alive=1,
            )
        else:
            print(f' |- it is a duplicate!')
            self.gadb.write(
                atoms,
                name=workdir,
                mag=0,
                eV=ene_eV,
                mated=0,
                done=1,
                alive=0,
            )
