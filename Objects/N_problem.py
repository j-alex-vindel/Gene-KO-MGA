from N_individual import N_Individual
import numpy as np
from Met_Net import MN
from typing import Type
from LP_MIBLP import Model_LP,Solve_LP

MetNet = Type[MN]

class LP_Problem:
    def __init__(self,metnet:MetNet=None,K:int=None,expand:bool=True):
        self.metnet = metnet
        self.K = K
        self.model = Model_LP(self.metnet)
        self.expand = expand

    def generate_individual(self):
        individual = N_Individual()
        individual.gene,individual.genesize,individual.wts = self.__genegenerator()
        return individual
    
    def calculate_objectives(self,individual):
        if self.expand: 
            individual.objectives = Solve_LP(network=self.metnet,model=self.model,yjs=individual.genbin)
        else:
            individual.objectives = Solve_LP(network=self.metnet,model=self.model,yjs=individual.genbin)


    def __genegenerator(self):
        flag = False
        while flag == False:
            if self.metnet.KO == None:
                gene = tuple(np.random.choice(self.metnet.M,self.K))
            else:
                gene = tuple(np.random.choice(self.metnet.KO,self.K))
            if len(set(gene)) == len(gene):
                flag=True
        return gene,len(self.metnet.M),(self.metnet.FBA[self.metnet.biomass],self.metnet.FBA[self.metnet.chemical])