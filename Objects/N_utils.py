from N_population import N_Pop
import random

class GA_utils:
    def __init__(self,problem,num_ind:int=30,num_par_tour:int=2,tour_prob:float=.9,mut_rate:float=1.0):
        self.num_ind = num_ind
        self.problem = problem
        self.num_par_tour = num_par_tour
        self.tour_prob = tour_prob
        self.mut_rate = mut_rate

    def create_init_population(self):
        population = N_Pop()
        for _ in range(self.num_ind):
            individual = self.problem.generate_individual()
            self.problem.calculate_objectives(individual)
            population.append(individual)
        return population
    
    def __crossover(self,individual1,individual2):
        child1 = self.problem.generate_individual()
        child2 = self.problem.generate_individual()
        A = list(individual1.gene)
        B = list(individual2.gene)
        k = random.randint(0,len(A)-1)
        a = A.pop(k)
        b = B.pop(k)
        A.extend([b])
        B.extend([a])
        child1.gene = tuple(A)
        child2.gene = tuple(B)
        return child1,child2 