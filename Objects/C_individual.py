from itertools import permutations


class Individual(object):
    def __init__(self):
        self.rank = None
        self.crowding_distance = None
        self.domination_count = None
        self.dominated_solutions = None
        self.gene = None
        self.objectives = None

    def __eq__(self,other):
        if isinstance(self,other.__class__):
            return self.gene == other.gene
        return False
    
    def dominates(self,other_individual):
        and_condition = True
        or_condition = False
        for first,second in zip(self.objectives,other_individual.objectives):
            and_condition = and_condition and first <= second
            or_condition = or_condition or first < second
        return (and_condition and or_condition)
    

class Ind(object):
    def __init__(self):
        self.gene:tuple= None
        self.objectives:tuple= None
        self.genesize:list = None
    
    @property
    def gene(self):
        return self._gene
    @gene.setter
    def gene(self,kindex:tuple):
        self._gene = kindex

    @property
    def genes(self):
        if self._gene == None:
            return None
        else:
            return list(permutations(self._gene))
    
    # Gene Size
    @property
    def genesize(self):
        return self._genesize
    @genesize.setter
    def genesize(self,value):
        self._genesize = value
    
    # Gene Binary Expresion
    @property
    def gene_bin(self):
        if (self._genesize == None) or (self._gene == None):
            return None
        else:
            return [0 if i in self._gene else 1 for i in range(self._genesize)]
        
    def __len__(self):
        return len(self._gene)