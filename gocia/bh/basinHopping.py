
from ase.db import connect
from gocia.interface import Interface
from ase.io import read, write

from gocia.bh import pass_Metropilis


class BasinHopping:
    def __init__(
        self,
        struct_curr = None,
        temperature = 300, #K
        substrate=None,
        subsPot=None,
        projname='bh',
        zLim=None,
        chemPotDict=None,
        convergeCrit=None,
        func_bias=None,
    ):

        self.temperature = temperature
        
        if type(struct_curr) is str:
            self.struct_curr = read(struct_curr)
        else:
            self.struct_curr = struct_curr

        if type(substrate) is str:
            self.substrate = read(substrate)
        else:
            self.substrate = substrate

        self.bhdb = connect(f'{projname}.db')

        if subsPot is not None:
            self.subsPot = subsPot
        else:
            try:
                self.subsPot = self.substrate.get_potential_energy()
            except:
                self.subsPot = 0

        if zLim is not None:
            self.zLim = zLim
        if zLim is None:
            self.zLim = Interface(self.substrate, self.substrate).zLim

        if chemPotDict is not None:
            self.chemPotDict = chemPotDict
            self.gc = True
        else:
            self.gc = False

        if convergeCrit is None:
            self.convergeCrit = 100
        else:
            self.convergeCrit = convergeCrit

        if func_bias:
            self.func_bias = func_bias
        else:
            self.func_bias = None

        self.rej_streak = 0
        self.gm_energy = self.struct_curr.get_potential_energy()
        self.gm_id = 0


    def __len__(self):
        return len(self.bhdb)
    
    def get_ID(self, condString):
        tmp = []
        for r in self.bhdb.select(condString):
            tmp.append(r.id)
        return tmp

    def get_valueOf(self, valueString, idList):
        tmp = []
        for i in idList:
            tmp.append(
                self.bhdb.get(id=i)[valueString]
            )
        return tmp
    
    def get_GMrow(self):
        if self.gc:
            eneList = self.get_valueOf('grandPot', self.get_ID('done=1'))
        else:
            eneList = self.get_valueOf('eV', self.get_ID('done=1'))
        self.Emin = min(eneList)
        row_min = None

        if self.gc:
            for r in self.bhdb.select(grandPot=self.Emin):
                row_min = r
        else:
            for r in self.bhdb.select(eV=self.Emin):
                row_min = r
        return row_min
    
    def is_converged(self):
        if len(self) == 0:
            return False
        
        gmid = self.get_GMrow().id
        if len(self) - gmid > self.convergeCrit:
            return True
        else:
            return False

    def calc_grandPot(self, atoms, epot = None):
    
        if epot is not None:
            myPot = epot - self.subsPot
        else:
            myPot = atoms.get_potential_energy() - self.subsPot
        adsSymbol = Interface(atoms, self.substrate).\
            get_adsAtoms().get_chemical_symbols()
        for s in adsSymbol:
            if s in self.chemPotDict:
                myPot -= self.chemPotDict[s]
                
        if self.func_bias:
            myPot += self.func_bias(atoms)
            
        return myPot
    
    def add_sample(self, struct_new, forced=False):

        if len(self) == 0:
            self.bhdb.write(
                struct_new,
                # mag=mag,
                eV = struct_new.get_potential_energy(),
                grandPot = self.calc_grandPot(struct_new),
                done=1,
                accept=1,
            )
            print(f'INIT    {struct_new.get_chemical_formula():12s}  {struct_new.get_potential_energy():6.3f}\n')
            return True

        if self.gc:
            delta_e =self.calc_grandPot(struct_new) - self.calc_grandPot(self.struct_curr)
        else:
            delta_e = struct_new.get_potential_energy() - self.struct_curr.get_potential_energy()

        # try:
        #     mag = self.struct_new.get_magnetic_moment()
        # except:
        #     mag = 0

        if pass_Metropilis(delta_e, self.temperature) or forced:
            self.struct_curr = struct_new
            self.bhdb.write(
                struct_new,
                # mag=mag,
                eV = struct_new.get_potential_energy(),
                grandPot = self.calc_grandPot(struct_new),
                done=1,
                accept=1,
            )
            self.rej_streak = 0
            if forced:
                print('FORCED', end='')
            else:
                print('ACCEPT', end='')
            print(f'  {struct_new.get_chemical_formula():12s}  {struct_new.get_potential_energy():6.3f}  |  CUR  {self.struct_curr.get_chemical_formula():12s}  {self.struct_curr.get_potential_energy():6.3f}\n')
            return True
        else:
            self.bhdb.write(
                struct_new,
                # mag=mag,
                eV = struct_new.get_potential_energy(),
                grandPot = self.calc_grandPot(struct_new),
                done=1,
                accept=0,
            )
            self.rej_streak += 1
            print(f'REJECT  {struct_new.get_chemical_formula():12s}  {struct_new.get_potential_energy():6.3f}  |  CUR  {self.struct_curr.get_chemical_formula():12s}  {self.struct_curr.get_potential_energy():6.3f}\n')
            return False
        



        
    
