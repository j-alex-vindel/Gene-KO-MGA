import sys 
import os
from typing import List,NewType,Type

import gurobipy as gp
from gurobipy import GRB
import copy
import random
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import pandas as pd
from N_individual import N_Individual


M_N = NewType('metnet',object)
IND = Type[N_Individual]
L = List[int]
FBA = NewType('Flux Vector',List[float])
WT = NewType('Wildtype',bool)
MT = NewType('Mutant',bool)
POP = NewType('POP',List[object])

def strainselector(strain:str=None) -> M_N:
    sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Metabolic_Networks'))) 
    
    match strain:
        case 'yeast':
            from YEAST import M_Yeast
            met = M_Yeast
        case 'ijo':
            from IJO import MN_ijo1366
            met = MN_ijo1366
        case 'iaf':
            from IAF import MN_iaf1260
            met = MN_iaf1260
        case 'ijr':
            from IJR import MN_ijr904
            met = MN_ijr904
        case 'ijrmomo':
            from IJRMOMO import MN_MOMO_ijr
            met = MN_MOMO_ijr
        case 'ijrmomod':
            from IJRMOMOD import MN_MOMO_D
            met = MN_MOMO_D
        case _:
            raise Exception("No Metabolic Network's been selected")
    return met

def plotmatrix(Matrix1=None,Matrix2=None,filename:str=None,title:str=None,filenames:list=None,save:bool=False):
     
     if filename != None:
        fig, axs = plt.subplots(1,2)
        ax1 = axs[0]
        ax2 = axs[1]

        ax1.spy(Matrix1,markersize=1,precision=.0001)
        ax1.set_xlabel('Space Explored')
        ax2.spy(Matrix2,markersize=1,precision=.0001,color='green')
        ax2.set_xlabel('Fitness')
        fig.suptitle(title)
        if save:
            plt.savefig(filename,format='jpg')
            filenames.append(filename)
            plt.close(fig)
        else:
            plt.show()
     else:
          pass


def make_gif(images:list=None,duration:int=250,filename:str=None):
     frames = [Image.open(image) for image in images]
     frame_one = frames[0]
     frame_one.save(filename,format='GIF',append_images=frames,save_all=True,duration=duration,loop=5)


def geneid(ind:IND=None,mnet:M_N=None):
     return [mnet.Rxn[i] for i in ind.gene]

def set_constructor(L) -> L:
    if L != None:
        return [_ for _ in range(len(L))]

def draw_fit_generation(y_list,title):
    x_list = np.arange(1, len(y_list)+1)  # create a numpy list from 1 to the numbers of generations

    plt.plot(x_list, y_list)

    plt.title(f"{title}")
    plt.xlabel("Generations")
    plt.ylabel("Fitness")

    plt.show()

def WT_FBA(mn,wt:WT=True,mt:MT=False) -> FBA:
    LB_wt = copy.deepcopy(mn.LB)
    UB_wt = copy.deepcopy(mn.UB)

    if wt and not mt:
        obj = mn.biomass
        FVA = False
    elif mt and not wt:
        obj = mn.chemical
        FVA=True
    else:
        raise Exception('Both wildtype and mutant cannot be TRUE or FALSE at the same time')
    
    fba = gp.Model()
    fba_v = fba.addVars(mn.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='vs')
    vs = [fba_v[i] for i in mn.M]

    fba.setObjective(1*fba_v[obj],GRB.MAXIMIZE)
    fba.addMConstr(mn.S,vs,'=',mn.b,name='Stoi')
    fba.addConstrs((LB_wt[j] <= fba_v[j] for j in mn.M), name='LBwt')
    fba.addConstrs((UB_wt[j] >= fba_v[j] for j in mn.M), name='UBwt')
    if FVA:
        fba.addConstr((fba_v[mn.biomass] >= mn.minprod), name='minprod')

    fba.Params.OptimalityTol = mn.infeas
    fba.Params.IntFeasTol = mn.infeas
    fba.Params.FeasibilityTol = mn.infeas
    fba.Params.OutputFlag = 0
    fba.optimize()

    if fba.status == GRB.OPTIMAL:
        fbavs =  [fba.getVarByName('vs[%s]'%a).x for a in mn.M] 
    else:
        fbavs = [-2000 for _ in mn.M]

    return fbavs

def save_population(population:POP=None,filename:str=None,num_gen:int=None,gen_size:int=None,mut_rate:float=None):
    res = {'Gene':[],
           'Bio':[],
           'Chem':[],
           'Fit':[],
           'NumGen':[],
           'MutRate':[],
           'GenSize':[]}
    
    for i in population:
        res['Gene'].append(i.gene)
        res['Bio'].append(i.objectives[0])
        res['Chem'].append(i.objectives[1])
        res['Fit'].append(i.fitness)
        res['NumGen'].append(num_gen)
        res['GenSize'].append(gen_size)
        res['MutRate'].append(mut_rate)

    df = pd.DataFrame.from_dict(res)

    df.to_csv(filename,index=False)