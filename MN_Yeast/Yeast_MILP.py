import sys
import os
from MetNet_YEAST import M_YEAST
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Objects')))

from MILP import MILP_MODEL,MILP_SOLVE

mn = M_YEAST

ys = [0 if i in [1,4,6] else 1 for i in mn.M]

milp = MILP_MODEL(network=mn)

solve = MILP_SOLVE(network=mn,model=milp,k=2,yjs=ys)

print(solve)