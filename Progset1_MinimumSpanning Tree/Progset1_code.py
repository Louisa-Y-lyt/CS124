import random
import numpy as np

class Graph:
    
    def __init__(self, vertices):
        self.V = vertices # No. of vertices
        self.graph = [] # to store graph
        self.parent = np.arange(0,self.V,1) 
        self.rank = [0] * self.V
        
    def set_graph(self, dimension):
        c = [3.15, 2.0, 1.82, 1.79]
        if dimension == 0:
            bar = 1
            if self.V > 5000:
                bar = c[0]*self.V**(-1)
            for i in range(self.V):
                res = np.random.rand(self.V-i-1)
                p = np.arange(i+1,self.V,1)
                for u, value in zip(p[res<bar],res[res<bar]):
                    self.graph.append([i, u, value])

        if dimension == 2:
            bar = 2
            if self.V > 5000:
                bar = c[1]*self.V**(-1/2)
            for i in range(self.V):
                x1 = np.random.rand(self.V-i-1)
                x2 = np.random.rand(self.V-i-1)
                res = (x1**2+x2**2)**0.5
                p = np.arange(i+1,self.V,1)
                for u, value in zip(p[res<bar],res[res<bar]):
                    self.graph.append([i, u, value])
 
        if dimension == 3:
            bar = 3
            if self.V > 5000:
                bar = c[2]*self.V**(-1/3)
            for i in range(self.V):
                x1 = np.random.rand(self.V-i-1)
                x2 = np.random.rand(self.V-i-1)
                x3 = np.random.rand(self.V-i-1)
                res = (x1**2+x2**2+x3**2)**0.5
                p = np.arange(i+1,self.V,1)
                for u, value in zip(p[res<bar],res[res<bar]):
                    self.graph.append([i, u, value])

        if dimension == 4:
            bar = 4
            if self.V > 5000:
                bar = c[3]*self.V**(-1/4)
            for i in range(self.V):
                x1 = np.random.rand(self.V-i-1)
                x2 = np.random.rand(self.V-i-1)
                x3 = np.random.rand(self.V-i-1)
                x4 = np.random.rand(self.V-i-1)
                res = (x1**2+x2**2+x3**2+x4**2)**0.5
                p = np.arange(i+1,self.V,1)
                for u, value in zip(p[res<bar],res[res<bar]):
                    self.graph.append([i, u, value])          
    
    def find(self, i):
        if self.parent[i] != i:
            self.parent[i] = self.find(self.parent[i])
        return self.parent[i]
    
    def link(self, i, j):
        if self.rank[i] > self.rank[j]:
            self.parent[j] = i
        if self.rank[i] < self.rank[j]:
            self.parent[i] = j
        if self.rank[i] == self.rank[j]:
            self.parent[j] = i
            self.rank[i] += 1
     
    def union(self,i,j):
        self.link(self.find(i), self.find(j))

    def kruskal(self):
        mst = []
        self.graph = sorted(self.graph, key=lambda item: item[2])
        for edge in self.graph:
            if self.find(edge[0]) != self.find(edge[1]):
                mst.append(edge[2])
                self.union(edge[0], edge[1])
        return np.sum(mst)




def runMST(flag = 0, numpoints = None, numtrials = None, dimension = None):
    # flag for running time record
    import random
    import numpy as np
    import time
    if dimension not in [0,2,3,4]:
        raise Exception("Wrong dimension parameter")
        return
    if flag != 0: 
        start = time.time()
    avg = []
    for t in range(numtrials):
        seed = t * dimension
        random.seed(seed)

        g = Graph(numpoints)
        g.set_graph(dimension)
        result = g.kruskal()
        avg.append(result)    
    average = np.mean(avg)
    if flag != 0: 
        end = time.time()
        print('The running time per trial is: {0} ms'.format((end-start) * 1000/numtrials)) 
    return [average, numpoints, numtrials, dimension]




if __name__ == "__main__":
    
    import sys

    print(runMST(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])))





