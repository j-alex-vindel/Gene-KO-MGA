from typing import List,NewType
import numpy as np
from Func import set_constructor,WT_FBA

S_M = List[List[int]]
L_B = List[int]
U_B = List[int]
RXN = List[str]
MET = List[str]
KOS = List[int]
B_M = NewType('Big M',int)
TGT = NewType('Target',float)
IDX = NewType('Index',int)


class MN:
    def __init__(self,
                 S:S_M=None,
                 LB:L_B=None,
                 UB:U_B=None,
                 Rxn:RXN=None,
                 Met:MET=None,
                 KO:KOS=None,
                 Name:str=None,
                 biomass:IDX=None,
                 chemical:IDX=None,
                 infeas:float=1e-6,
                 time_limit:int=1000,
                 BM:B_M=1000,
                 ):
        self.S = S
        self.LB = LB
        self.UB = UB
        self.S = S
        self.LB = LB
        self.UB = UB
        self.Rxn = Rxn
        self.Met = Met
        self.KO = KO
        self.Name = Name
        self.biomass = biomass
        self.chemical = chemical
        self.infeas = infeas
        self.time_limit = time_limit
        self.BM = BM
        self.M = set_constructor(self.Rxn)
        self.N = set_constructor(self.Met)
        self.b = np.array([0 for i in self.N])
        self.c = np.array([1 if i == self.biomass else 0 for i in self.M])
        self.tgt = .5
        self.FBA = WT_FBA(self)
        self.FVA = WT_FBA(self,wt=False,mt=True)


    @property
    def tgt(self):
        return self._tgt
    @tgt.setter
    def tgt(self,tgt:TGT=.5):
        self._minprod = None
        self._tgt = tgt
    @property
    def minprod(self):
        if self._minprod is None:
            self._minprod = self._tgt*self.FBA[self.biomass]
        return self._minprod