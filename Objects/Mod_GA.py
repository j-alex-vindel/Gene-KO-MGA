from N_individual import N_Individual
from Met_Net import MN
from LP_MIBLP import Model_LP, Solve_LP
import statistics
import numpy as np
import random
from typing import List,Type,Tuple,NewType
import copy
from tqdm import tqdm


IND = Type[N_Individual]
MetNet = Type[MN]
Model = Type[Model_LP]
NumGen = NewType('NumGen',int)
MutRate = NewType('MutRate',float)

def create_children(population,numpar=2,prob=.4):
     children = []
     while len(children) < len(population):
          parent1 = tournament(population,num_par=numpar,prob=prob)
          parent2 = parent1
          while parent1 == parent2:
               parent2 = tournament(population,num_par=numpar,prob=prob)
          child1, child2 = crossover(parent1,parent2)
          
          
          children.append(child1)
          children.append(child2)
     
     
     
     return children


def tournament(population,num_par,prob):
     participants = random.sample(population,num_par)
     best = None
     for participant in participants:
          if best is None or ((participant.cost >=0 )and (choose_w_prob(prob))):
               best = participant
     return best


def genegenerator(metnet=None,K=None):
        flag = False
        while flag == False:
            if metnet.KO == None:
                gene = tuple(np.random.choice(metnet.M,K))
            else:
                gene = tuple(np.random.choice(metnet.KO,K))
            if len(set(gene)) == len(gene):
                flag=True
        return gene, len(metnet.M),(metnet.FBA[metnet.biomass],metnet.FBA[metnet.chemical])


def crossover(parent1:IND,parent2:IND):
    child1 = N_Individual()
    child2 = N_Individual()
    A = copy.deepcopy(list(parent1.gene))
    B = copy.deepcopy(list(parent2.gene))
    k = random.randint(0,len(A)-1)
    a = A.pop(k)
    b = B.pop(k)
    A.extend([b])
    B.extend([a])
    child1.gene = tuple(A)
    child2.gene = tuple(B)
    child1 = copy.deepcopy(child1)
    child2 = copy.deepcopy(child2)
    child1.genesize,child1.wts = parent1.genesize,parent1.wts
    child2.genesize,child2.wts = parent2.genesize,parent2.wts
    return child1,child2

def choose_w_prob(prob:float=None):
     if random.random() >= prob:
          return True
     return False

def mutate(child:IND=None,mut_rate:float=.7):
     
     chance = random.random()
     L = list(child.gene)
     if (chance <= mut_rate) or (len(set(L)) != len(L)):
           ri = random.randint(0,len(L)-1)
           flag = False
           while flag == False:
                new_value = random.randint(0,child.genesize-1)
                if new_value not in L:
                     flag = True
                     L[ri] = new_value
     child.gene = tuple(L)

def mutate_2(child:IND=None,mut_rate:float=.7,mut_scope:int=5):
     chance = random.random()
     L = list(child.gene)
     if (chance <= mut_rate) or (len(set(L)) != len(L)):
          t1,t2 = L[0] + mut_scope, L[1] +mut_scope
          l1,l2 = L[0] - mut_scope, L[1] - mut_scope
          if t1 >child.genesize:
               t1 = child.genesize
          if t2 > child.genesize:
               t2 = child.genesize
          if l1 < 0:
               l1 =0
          if l2 <0:
               l2 =0
          flag = False
          while flag ==False:
               v1 = random.randint(l1,t1)
               v2 = random.randint(l2,t2)
               NL = [v1,v2]
               if NL != L:
                    flag=True
                    L = NL
     child.gene = tuple(L)
     
def evolve(pop_size:int=None,K:int=None,metnet:MetNet=None,model:Model=None,num_gen:NumGen=None,numpar:int=None,crossprob:float=None,mutrate:MutRate=.7):
     B_Matrix = np.zeros([len(metnet.M),len(metnet.M)])
     C_Matrix = np.zeros([len(metnet.M),len(metnet.M)])
     F_Matrix = np.zeros([len(metnet.M),len(metnet.M)])
     Explored = np.zeros([len(metnet.M),len(metnet.M)])
     
     ppop = [N_Individual() for _ in range(pop_size)]
     ppop_full = []
     AVGs = []
     HIGH_gen = []
     for individual in ppop:
          individual.gene,individual.genesize,individual.wts = genegenerator(metnet=metnet,K=K) 
          individual.objectives = Solve_LP(network=metnet,model=model,yjs=individual.genebin)
          for gene in individual.genes:
               if Explored[gene] == 0:
                    B_Matrix[gene] = individual.objectives[0]
                    C_Matrix[gene] = individual.objectives[1]
                    F_Matrix[gene] = individual.fitness
                    Explored[gene] = 1
     
     avg = statistics.fmean([x.fitness for x in ppop])
     high = max([x.fitness for x in ppop])
     AVGs.append(avg)
     HIGH_gen.append(high)
     for it in tqdm(range(num_gen)):
          cpop = create_children(population=ppop,numpar=numpar,prob=crossprob)
          for child in cpop:
               mutate(child=child,mut_rate=mutrate)
               if Explored[child.gene] == 0:
                  child.objectives = Solve_LP(network=metnet,model=model,yjs=child.genebin)
                  for gene in child.genes:
                       if Explored[gene] == 0:
                            B_Matrix[gene] = child.objectives[0]
                            C_Matrix[gene] = child.objectives[1]
                            F_Matrix[gene] = child.fitness
                            Explored[gene] = 1
               else:
                  child.objectives = B_Matrix[child.gene],C_Matrix[child.gene]
          ppop += cpop
          ppop_full += ppop
          ppop = sorted(ppop,key=lambda x:x.fitness,reverse=True)

          ppop = ppop[0:pop_size]
          avg = statistics.fmean([x.fitness for x in ppop])
          AVGs.append(avg) 
          high = max([x.fitness for x in ppop])
          HIGH_gen.append(high)
     ppop_full = sorted(ppop_full,key=lambda x:x.fitness, reverse=True)
     
     return ppop, ppop_full, AVGs,HIGH_gen

