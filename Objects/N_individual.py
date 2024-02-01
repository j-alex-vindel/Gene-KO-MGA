from itertools import permutations

class N_Individual(object):
    def __init__(self):
        self.gene:tuple = None
        self.objectives:tuple = None
        self.genesize:int= None 
        self.wts:tuple=None
        self.genesid:list=None
    def __len__(self):
        return len(self._gene)
    # Properties and its setters ###############
    @property
    def gene(self):
        return self._gene
    @gene.setter
    def gene(self,value:tuple):
        self._gene = value 
        
    @property
    def genesize(self):
        return self._genesize
    @genesize.setter
    def genesize(self,value:int):
        self._genesize = value 
    
    @property
    def wts(self):
        return self._wts 
    @wts.setter
    def wts(self,value):
        self._wts = value 
    
    @property
    def genesid(self):
        if self._gene == None:
            return None
        else:
            return self._genesid
    @genesid.setter
    def genesid(self,value):
        self._genesid = value
    #  Just Properties #####################
    @property
    def genebin(self):
        if (self._genesize == None) or (self._gene == None):
            return None
        else:
            return [0 if i in self._gene else 1 for i in range(self._genesize)]
        
    @property
    def genes(self):
        if self._gene == None:
            return None
        else:
            return list(permutations(self._gene))
    
    @property
    def cost(self):
        if (self.objectives == None) or (self._wts == None):
            return None
        else:
            return sum([1000*(self.objectives[0]),(self.objectives[1])])
        
    @property
    def fitness(self):
        if self.objectives == None:
            return None
        else:
            f1,f2 = self.objectives[0],self.objectives[1]
            if f1 ==0:
                v1 = 0
            else:
                v1 = 1/f1
            if f2 == 0:
                v2 = 0
            else:
                v2 = 1/f2
            return v1+v2
  