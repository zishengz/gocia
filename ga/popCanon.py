
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
        substrate = None,
        gadb = None,
        popSize = 20,
        zLim = None,
        compParam = None,
        matingMethod = None,
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

        if chemPotDict is not None:
            self.chemPotDict = chemPotDict

        self.iniSize = len(self)

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
            aliveList = [x for _,x in sorted(zip(fitnessList, aliveList))]
            deadList = aliveList[:-self.popSize]
            for d in deadList:
                print('REST IN PEACE, %i!'%d)
                self.gadb.update(d, alive=0)

    def is_uniqueInPop(self, atoms):
        '''
        Check similarity against the current population
        '''
        aliveList = self.get_ID('alive=1')
        aliveList = [self.gadb.get(id=n).toatoms() for n in aliveList]
        isUnique = True
        for a in aliveList:
            if srtDist_similar_zz(atoms, a):
                isUnique = False
                break
        return isUnique

    def gen_offspring(self, mutRate=0.4):
        kid = None
        mater, pater = 0, 0
        while kid is None:
            mater, pater = self.choose_parents()
            a1 = self.gadb.get(id=mater).toatoms()
            a2 = self.gadb.get(id=pater).toatoms()
            surf1 = Interface(a1, self.substrate, zLim=self.zLim)
            surf2 = Interface(a2, self.substrate, zLim=self.zLim)
            kid = crossover_snsSurf_2d(surf1, surf2, tolerance=0.75)
        print('PARENTS: %i and %i'%(mater, pater))
        if srtDist_similar_zz(a1, a2)\
            or srtDist_similar_zz(a1, kid.get_allAtoms())\
            or srtDist_similar_zz(a2, kid.get_allAtoms()):
            print(' |- TOO SIMILAR!')
            mutRate = 1
        if np.random.rand() < mutRate:
            print(' |- MUTATION!')
            kid.transMut()
            kid.rattleMut()
        self.gadb.update(mater, mated=self.gadb.get(id=mater).mated+1)
        self.gadb.update(pater, mated=self.gadb.get(id=pater).mated+1)
        return kid

    def add_vaspResult(self, vaspdir='.'):
        import os
        cwdFiles = os.listdir(vaspdir)
        if 'OSZICAR' in cwdFiles and 'BADSTRUCTURE' not in cwdFiles:
            info = [l for l in open('%s/OSZICAR'%vaspdir).readlines() if 'F' in l]
            if len(info) > 0:
                info = info[-1].split()
                mag = eval(info[-1])
                ene_eV = eval(info[4])
                s = read('%s/CONTCAR'%vaspdir)
                s.wrap()
                print('\nA CHILD IS BORN with G = %.3f eV'%(ene_eV))
                if self.is_uniqueInPop(s):
                    if ene_eV < self.get_GMrow()['ene_eV']:
                        print(' |- it is the new GM!')
                    self.gadb.write(
                        s,
                        mag     = mag,
                        eV      = ene_eV,
                        mated   = 0,
                        done    = 1,
                        alive   = 1
                    )
                else:
                    print(' |- it is a duplicate!')
                    self.gadb.write(
                        s,
                        mag     = mag,
                        eV      = ene_eV,
                        mated   = 0,
                        done    = 1,
                        alive   = 0
                    )





    