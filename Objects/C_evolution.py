from C_utils import NSGA_UTILS
from C_population import Population
from tqdm import tqdm

class Evolution:
    def __init__(self,problem, num_gen=100, num_ind=30, num_par_tour=2,tour_prob=.9,mutation_rate=1):
        self.utils = NSGA_UTILS(problem=problem, num_ind=num_ind,num_par_tour=num_par_tour,tour_prob=tour_prob,mutation_rate=mutation_rate)
        self.population = None
        self.num_gen = num_gen
        self.on_gen_finished = []
        self.num_ind = num_ind

    def evolve(self):
        self.population = self.utils.create_initial_population()
        self.utils.fast_nondominated_sort(self.population)
        for front in self.population.fronts:
            self.utils.calculate_crowding_distance(front)
        children = self.utils.create_children(self.population)
        returned_population = None
        for i in tqdm(range(self.num_gen)):
            self.population.extend(children)
            self.utils.fast_nondominated_sort(self.population)
            new_population = Population()
            front_num = 0
            while len(new_population) + len(self.population.fronts[front_num]) <= self.num_ind:
                self.utils.calculate_crowding_distance(self.population.fronts[front_num])
                new_population.extend(self.population.fronts[front_num])
                front_num += 1
            self.utils.calculate_crowding_distance(self.population.fronts[front_num])
            self.population.fronts[front_num].sort(key=lambda individual: individual.crowding_distance, reverse=True)
            new_population.extend(self.population.fronts[front_num][0:self.num_ind - len(new_population)])
            returned_population = self.population
            self.population = new_population
            self.utils.fast_nondominated_sort(self.population)
            for front in self.population.fronts:
                self.utils.calculate_crowding_distance(front)
            children = self.utils.create_children(self.population)
        return returned_population.fronts 