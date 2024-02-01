from N_individual import N_Individual
import numpy as np
import random
from Func import strainselector, plotmatrix,make_gif,geneid
from LP_MIBLP import Model_LP, Solve_LP
import copy
from typing import Type
import matplotlib.pyplot as plt
import dill

import json

IND = Type[N_Individual]

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
     

#####m=--------------------------====####################3
version = 'A4'
maxit = 10
Npop = 10
mn = strainselector('ijo')
BIO_MATRIX = np.zeros([len(mn.M),len(mn.M)])
CHE_MATRIX = np.zeros([len(mn.M),len(mn.M)]) 
FIT_MATRIX = np.zeros([len(mn.M),len(mn.M)])
FIT_POS_MATRIX = np.zeros([len(mn.M),len(mn.M)])
mut_rate = .9
model = Model_LP(network=mn)
graphs = False
POP = []
ppop = [N_Individual() for _ in range(Npop)]

matrix = {}
filenames = []
tot_individuals = 0

if graphs:
     title = f"Generation 0 \n Individuals {tot_individuals} \n Mutation Rate:{int(mut_rate*100)}%"
     filename = f"../Figs/Gen_{version}0.jpg"
     plotmatrix(Matrix1=FIT_MATRIX,Matrix2=FIT_POS_MATRIX,filename=filename,title=title,filenames=filenames)

for individual in ppop:
    individual.gene,individual.genesize,individual.wts = genegenerator(metnet=mn,K=2)
    individual.objectives = Solve_LP(network=mn,model=model,yjs=individual.genebin)
    matrix[individual.gene] = individual.objectives
    
    for gene in individual.genes:
        if BIO_MATRIX[gene] == 0:
            BIO_MATRIX[gene] = individual.objectives[0]
            CHE_MATRIX[gene] = individual.objectives[1] 
            FIT_MATRIX[gene] = individual.cost
            if individual.cost >=0:
                 FIT_POS_MATRIX[gene] = individual.cost

tot_individuals += len(ppop)

if graphs:
     title = f"Generation 1 \n Individuals {tot_individuals} \n Mutation Rate:{int(mut_rate*100)}%"
     filename = f"../Figs/Gen_{version}1.jpg"
     plotmatrix(Matrix1=FIT_MATRIX,Matrix2=FIT_POS_MATRIX,filename=filename,title=title,filenames=filenames)

# print(np.transpose((BIO_MATRIX<0).nonzero()))
# print('Axis = 0')
# print(np.where(BIO_MATRIX.any(axis=(0,1))[0]))

# Matrix = np.empty([3,3],dtype=[('biomass','f4'),('chemical','f4')])

tot_individuals = 0
for it in range(maxit): 
    tot_individuals += len(ppop)
    cpop = create_children(population=ppop,numpar=2,prob=.4)


    print(f"\n{len(cpop)} New Individuals - Generation {it+1}")
    for child in cpop:
     
        mutate_2(child=child,mut_rate=mut_rate,mut_scope=6)
        if child.gene not in matrix:
            child.objectives = Solve_LP(network=mn,model=model,yjs=child.genebin)
            matrix[child.gene] = child.objectives
            for gene in child.genes:
                if BIO_MATRIX[gene] == 0:
                    BIO_MATRIX[gene] = child.objectives[0]
                    CHE_MATRIX[gene] = child.objectives[1]
                    FIT_MATRIX[gene] = child.cost
                    if individual.cost >=0:
                        FIT_POS_MATRIX[gene] = individual.cost
        else:
             child.objectives = matrix[child.gene]
    if graphs:
          title = f"Generation {it+2} \n Individuals {tot_individuals}\n Mutation Rate:{int(mut_rate*100)}%"
          filename = f"../Figs/Gen_{version}{it+2}.jpg"
          plotmatrix(Matrix1=FIT_MATRIX,Matrix2=FIT_POS_MATRIX,filename=filename,title=title,filenames=filenames)
    
    ppop += cpop
    POP += ppop
    ppop = sorted(ppop,key=lambda x:x.cost, reverse=True)
    ppop = ppop[0:Npop]


# print(np.transpose((BIO_MATRIX<0).nonzero()))
for individual in ppop:
     individual.genesid = geneid(ind=individual,mnet=mn)
     print(f"{individual.gene} -> {individual.objectives} -> {individual.cost} -> {individual.fitness} -> {individual.genesid}")


if graphs:
     gifname = f"../Figs/GIF_{version}.gif"
     make_gif(filename=gifname,images=filenames,duration=300)


# dillname = f"../Results/Population_GA{version}.pkl"
# with open(dillname,'wb') as f:
#      dill.dump(ppop,f)

print(len(matrix))
print(len(POP))
