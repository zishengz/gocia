
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

        self.popSize = popSize

        eneList = self.get_valueOf('eV', self.get_ID('done=1'))
        self.Emin = min(eneList)

    def initializeDB(self):
        for i in range(len(self.gadb)):
            self.gadb.update(i+1, mated=0, alive=1)

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
            aliveList = [x for _,x in sorted(zip(fitnessList, aliveList))]
            deadList = aliveList[:-self.popSize]
            for d in deadList:
                print('REST IN PEACE, %i!'%d)
                self.gadb.update(d, alive=0)

    def gen_offspring(self, mutRate=0.3):
        kid = None
        mater, pater = 0, 0
        while kid is None:
            mater, pater = self.choose_parents()
            a1 = self.gadb.get(id=mater).toatoms()
            a2 = self.gadb.get(id=pater).toatoms()
            surf1 = Interface(a1, self.substrate, zLim=self.zLim)
            surf2 = Interface(a2, self.substrate, zLim=self.zLim)
            kid = crossover_snsSurf_2d(surf1, surf2, tolerance=0.8)
        print('PARENTS: %i and %i'%(mater, pater))
        if srtDist_similar_zz(a1, a2):
            print(' |- TOO SIMILAR!')
            mutRate *= 2
        if np.random.rand() < mutRate:
            print(' |- MUTATION!')
            kid.rattleMut()
        self.gadb.update(mater, mated=self.gadb.get(id=mater).mated+1)
        self.gadb.update(pater, mated=self.gadb.get(id=pater).mated+1)
        return kid

    def add_vaspResult(self, vaspdir='.'):
        import os
        if 'OSZICAR' in os.listdir(vaspdir):
            info = [l for l in open('%s/OSZICAR'%vaspdir).readlines() if 'F' in l]
            if len(info) > 0:
                info = info[-1].split()
                mag = eval(info[-1])
                ene_eV = eval(info[4])
                s = read('%s/CONTCAR'%vaspdir)
                s.wrap()
                print('\nA CHILD IS BORN with E = %.3f eV'%(ene_eV))
                self.gadb.write(
                    s,
                    mag     = mag,
                    eV      = ene_eV,
                    mated   = 0,
                    done    = 1,
                    alive   = 1
                )





    