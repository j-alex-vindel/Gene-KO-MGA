import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Objects')))
from C_evolution import Evolution
from C_problem import BILP
from Func import strainselector,gene2name

metnet = strainselector('yeast')

print(metnet.Name)
metnet.tgt = .5
problem = BILP(metnet=metnet,K=2)
evo = Evolution(problem=problem,num_gen=10,num_ind=5)

nga = evo.evolve()

func = [ i for i in nga]

obs = [gene2name(network=metnet,gene=i.gene) for i in func[0]]

print(obs)