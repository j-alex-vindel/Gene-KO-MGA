from C_population import Population
import random


class NSGA_UTILS:

    def __init__(self, problem,num_ind:int=30,num_par_tour:int=2,tour_prob:float=.9,mutation_rate:float=1.0):
        self.num_ind = num_ind
        self.problem = problem
        self.num_par_tour = num_par_tour
        self.tour_prob = tour_prob
        self.mutation_rate = mutation_rate
        

    def create_initial_population(self):
        population = Population()
        for _ in range(self.num_ind):
            individual = self.problem.generate_individual()
            self.problem.calculate_objectives(individual)
            population.append(individual)
        return population
    
    def fast_nondominated_sort(self, population):
        population.fronts = [[]]
        for individual in population:
            individual.domination_count = 0
            individual.dominated_solutions = []
            for other_individual in population:
                if individual.dominates(other_individual):
                    individual.dominated_solutions.append(other_individual)
                elif other_individual.dominates(individual):
                    individual.domination_count += 1
            if individual.domination_count == 0:
                individual.rank = 0
                population.fronts[0].append(individual)
        i = 0
        while len(population.fronts[i]) > 0:
            temp = []
            for individual in population.fronts[i]:
                for other_individual in individual.dominated_solutions:
                    other_individual.domination_count -= 1
                    if other_individual.domination_count == 0:
                        other_individual.rank = i + 1
                        temp.append(other_individual)
            i = i + 1
            population.fronts.append(temp)

    def calculate_crowding_distance(self, front):
        if len(front) > 0:
            solutions_num = len(front)
            for individual in front:
                individual.crowding_distance = 0

            for m in range(len(front[0].objectives)):
                front.sort(key=lambda individual: individual.objectives[m])
                front[0].crowding_distance = 10 ** 9
                front[solutions_num - 1].crowding_distance = 10 ** 9
                m_values = [individual.objectives[m] for individual in front]
                scale = max(m_values) - min(m_values)
                if scale == 0: scale = 1
                for i in range(1, solutions_num - 1):
                    front[i].crowding_distance += (front[i + 1].objectives[m] - front[i - 1].objectives[m]) / scale

    def crowding_operator(self, individual, other_individual):
        if (individual.rank < other_individual.rank) or \
                ((individual.rank == other_individual.rank) and (
                        individual.crowding_distance > other_individual.crowding_distance)):
            return 1
        else:
            return -1

    def create_children(self, population):
        children = []
        while len(children) < len(population):
            parent1 = self.__tournament(population)
            parent2 = parent1
            while parent1 == parent2:
                parent2 = self.__tournament(population)
            child1, child2 = self.__crossover(parent1, parent2)
            self.__mutate(child1)
            self.__mutate(child2)
            self.problem.calculate_objectives(child1)
            self.problem.calculate_objectives(child2)
            children.append(child1)
            children.append(child2)

        return children

    def __crossover(self, individual1, individual2):
        child1 = self.problem.generate_individual()
        child2 = self.problem.generate_individual()
        length = len(individual1.gene)
        ua = [random.randint(0,1) for _ in range(length)]
        for i in range(length):
            child1.gene[i] = ua[i]*individual1.gene[i] + (1-ua[i]*individual2.gene[i])
            child2.gene[i] = ua[i]*individual2.gene[i] + (1-ua[i]*individual1.gene[i])

        return child1, child2


    def __mutate(self, child):
        mutation_index = [random.choice(self.problem.metnet.M) for _ in range(self.mutation_rate)]
        for i, gene in enumerate(child.gene):
            if i not in mutation_index:
                child.gene[i] = gene
            else:
                child.gene[i] = 1-gene

  

    def __tournament(self, population):
        participants = random.sample(population.population, self.num_par_tour)
        best = None
        for participant in participants:
            if best is None or (
                    self.crowding_operator(participant, best) == 1 and self.__choose_with_prob(self.tour_prob)):
                best = participant

        return best

    def __choose_with_prob(self, prob):
        if random.random() <= prob:
            return True
        return False