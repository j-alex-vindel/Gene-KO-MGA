import gurobipy as gp
from gurobipy import GRB
from typing import List, Tuple,Type
from collections import namedtuple
from Met_Net import MN

import copy


Vector = List[float]
Model = Type[object]
Y = List[int]
M_N = Type[MN]
Result = namedtuple('Result',['MetNet','Strategy','Vs','Time','Soltype'])
RMILP = namedtuple('RMILP',['Cost','Chem','Biom'])
K = Type[int]


def Model_LP(network:M_N=None) -> Model:
    '''
    Returns a gurobi LP Model
    '''
    m = gp.Model('LP_Model')
    # Variables
    v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
    vs = [v[i] for i in network.M]
    # Objective
    m.setObjective((network.BM*v[network.biomass] + v[network.chemical]),GRB.MAXIMIZE)

    # Constraints

    m.addMConstr(network.S,vs,'=',network.b,name='S')

    m.update()

    model = m.copy()

    return model

def Solve_LP(network:M_N=None,model:Model=None,yjs:Y=None):
    '''
    Solves the LP based on a binary vector based on the ko strategy
    '''
    # retrieves variables
    v = [model.getVarByName('v[%s]'%a) for a in network.M]
    # Bounds
    model.addConstrs((network.LB[j]*yjs[j] <= v[j] for j in network.M), name='LB')
    model.addConstrs((v[j] <= network.UB[j]*yjs[j] for j in network.M), name='UB')
    # Add Parameters
    model.Params.OptimalityTol = network.infeas
    model.Params.IntFeasTol = network.infeas
    model.Params.FeasibilityTol = network.infeas
    model.Params.NodefileStart = 0.5
    # model.Params.Presolve = 0
    model.Params.OutputFlag = 0
    model.update()
    model.optimize()
     # Update
    # Solve
    if model.status == GRB.OPTIMAL:
        vs = [model.getVarByName('v[%d]'%j).x for j in network.M]

    elif model.status == GRB.TIME_LIMIT:
        vs = [model.getVarByName('v[%d]'%j).x for j in network.M]
    
    if model.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
        # print('Model status: *** INFEASIBLE or UNBOUNDED ***')
        vs = [-2000 for i in network.M]
    
    chem = vs[network.chemical]
    biom = vs[network.biomass]
    
    return  (biom,chem)


def MILP_MODEL(network:M_N=None) -> Model:
    '''
    Returns a gurobi model ready to be solved
    '''
    m = gp.Model("Single_Level_Reformulation")
    #Variables
    v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
    ys = m.addVars(network.M,lb=0,ub=1,vtype=GRB.INTEGER,name='ys')
    vs = [v[i] for i in network.M]
    # Dual Variables
    l = m.addVars(network.N,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='l')
    a1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='a1')
    b1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='b1')
    a2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='a2')
    b2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='b2')
    a = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='a')
    b = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='b')
    
    # Objective
    m.setObjective((1*v[network.chemical]),GRB.MAXIMIZE)
    # Constraints
    # Stoichimetric Constrs
    m.addMConstr(network.S,vs,'=',network.b,name='S')
    
    # Dual Objective
    m.addConstr((v[network.biomass] >= (sum(a1[j]*network.UB[j] - b1[j]*network.LB[j] for j in network.M)
     + sum(a2[j]*network.UB[j] - b2[j]*network.LB[j] for j in network.M))),name='DO')
    
    # Dual Constraints
    m.addConstrs((gp.quicksum(network.S.transpose()[i,j]*l[j] for j in network.N)
              - b[i]
              + a[i] - b2[i] + a2[i]
               == network.c[i] for i in network.M)
             ,name='SD')

    m.addConstr(v[network.biomass] >= network.minprod, name='min_growth')

    m.update()

    model = m.copy()

    return model

def MILP_SOLVE(network:M_N,model:Model=None,yjs:Y=None,k:K=None):
    '''
    Takes a y vector to solve the milp model
    '''

    # Retrieve variables from model
    v = [model.getVarByName('v[%s]'%a) for a in network.M]
    ys = [model.getVarByName('ys[%s]'%a) for a in network.M]
    l = [model.getVarByName('l[%s]'%a) for a in network.N]
    a = [model.getVarByName('a[%s]'%a) for a in network.M]
    a1 = [model.getVarByName('a1[%s]'%a) for a in network.M]
    a2 = [model.getVarByName('a2[%s]'%a) for a in network.M]
    b = [model.getVarByName('b[%s]'%a) for a in network.M]
    b1 = [model.getVarByName('b1[%s]'%a) for a in network.M]
    b2 = [model.getVarByName('b2[%s]'%a) for a in network.M]
    
    model.addConstrs((ys[i]==yjs[i] for i in network.M))
    # Add Linerarization Constraints
    # Linearization
    model.addConstrs((a1[j] <= network.BM*ys[j] for j in network.M),name='l1_a1')

    model.addConstrs((a1[j] >= - network.BM*ys[j] for j in network.M),name='l2_a1')

    model.addConstrs((a1[j] <= a[j] + network.BM*(1-ys[j]) for j in network.M),name='l3_a1')

    model.addConstrs((a1[j] >= a[j] - network.BM*(1-ys[j]) for j in network.M),name='l4_a1')

    model.addConstrs((b1[j] <= network.BM*ys[j] for j in network.M),name='l1_b1')

    model.addConstrs((b1[j] >= -network.BM*ys[j] for j in network.M),name='l2_b1')

    model.addConstrs((b1[j] <= b[j] + network.BM*(1-ys[j]) for j in network.M),name='l3_b1')

    model.addConstrs((b1[j] >= b[j] - network.BM*(1-ys[j]) for j in network.M),name='l4_b1')

    # Bounds
    model.addConstrs((network.LB[j]*ys[j] <= v[j] for j in network.M), name='LB')
    model.addConstrs((v[j] <= network.UB[j]*ys[j] for j in network.M), name='UB')

    model.addConstrs((network.LB[j] <= v[j] for j in network.M),name='lb')
    model.addConstrs((v[j] <= network.UB[j] for j in network.M),name='ub')

    # Knapsack
    if network.KO is not None:

        model.addConstr(sum(1-ys[j] for j in network.KO) == k,name='ksack')
        model.addConstrs((ys[j] == 1 for j in network.M if j not in network.KO),name='essen')
    elif network.KO is None:
        model.addConstr(sum(1-ys[j] for j in network.M) == k, name='ksackall')

    # Add Parameters
    model.Params.OptimalityTol = network.infeas
    model.Params.IntFeasTol = network.infeas
    model.Params.FeasibilityTol = network.infeas
    model.Params.NodefileStart = 0.5
    # model.Params.Presolve = 0
    model.Params.OutputFlag = 0
    model.update()
    model.optimize()

    # Update
    # Solve
    if model.status == GRB.OPTIMAL:
        chem = model.getObjective().getValue()
        vs = [model.getVarByName('v[%d]'%j).x for j in network.M]

    elif model.status == GRB.TIME_LIMIT:
        vs = [model.getVarByName('v[%d]'%j).x for j in network.M]
    
    if model.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
        # print('Model status: *** INFEASIBLE or UNBOUNDED ***')
        vs = [-2000 for i in network.M]
    
    chem = vs[network.chemical]
    biom = vs[network.biomass]
    
    return  [-biom,-chem]