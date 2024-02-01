from Mod_GA import evolve
from Func import strainselector,geneid,draw_fit_generation,save_population,plotmatrix
from LP_MIBLP import Model_LP

## ---- Input Parameters -----
num_gen = 15
pop_size = 10
strain = 'ijo'
mut_rate=.9
K=2
c_prob =.4
version = 'D'
outfilename = f"../Results/GA_Version{version}.csv"
## ---- Model Setup ----
Mn = strainselector(strain=strain)
Model = Model_LP(network=Mn)

## ---- Evolve ----
Pop,Pop_full,avgs,highs = evolve(pop_size=pop_size,K=K,metnet=Mn,model=Model,num_gen=num_gen,numpar=2,crossprob=c_prob,mutrate=mut_rate)

for i in Pop:
    print(f"{i.gene} -> {i.objectives} -> {i.fitness} -> {geneid(ind=i,mnet=Mn)}")
print(f" \n")



draw_fit_generation(avgs,'Avg Fitness through generations')

draw_fit_generation(highs,'Max Ftiness through generations')
# plotmatrix(Matrix1=explored,Matrix2=fit,title="Matrixs",filename='Matrix_G15_P10')

# save_population(population=Pop,filename=outfilename,num_gen=num_gen,gen_size=pop_size,mut_rate=mut_rate)