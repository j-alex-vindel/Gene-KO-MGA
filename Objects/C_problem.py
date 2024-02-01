from C_individual import Individual
import random
from LP_MIBLP import MILP_SOLVE,MILP_MODEL
from Met_Net import MN
from typing import Type
import numpy as np

M_N = Type[MN]
class BILP:
    def __init__(self,metnet:M_N=None,K:int=None,expand:bool=True):
        self.metnet = metnet
        self.K = K
        self.model = MILP_MODEL(self.metnet) 
        self.expand = expand

    def generate_individual(self):
        individual = Individual()
        individual.gene = list(np.random.choice([0,1],size=len(self.metnet.M)))
        return individual
    
    def calculate_objectives(self,individual):
        if self.expand:
            individual.objectives = MILP_SOLVE(network=self.metnet,model=self.model,k=self.K,yjs=individual.gene)
        else:
            individual.objectives = MILP_SOLVE(network=self.metnet,model=self.model,k=self.K,yjs=individual.gene)