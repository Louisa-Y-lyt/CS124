#!/usr/bin/env python
# coding: utf-8

# In[1]:


def read_file(input_file):
    with open(input_file, 'r') as f:
        rawdata = np.genfromtxt(f, dtype=int, delimiter=',')
    return rawdata


# In[8]:


def KK(input_seq):
    heap = []
    for num in input_seq:
        heapq.heappush(heap, -num)
    while len(heap) >= 2:
        a = heapq.heappop(heap)
        b = heapq.heappop(heap)
        heapq.heappush(heap, a-b)
        
    resd = -heapq.heappop(heap)
    return resd


# In[63]:


def random_solutionGen(n):
    return [random.randint(0, 1) * 2 - 1 for _ in range(n)]

def residue(input_seq, solution):
    return abs(np.dot(input_seq, solution))

def random_neighbors(solution):
    s = copy.deepcopy(solution)
    n = len(solution)
    i, j = random.sample(range(n), 2)
    s[j] = -s[j]
    if random.random() < 0.5:
        s[i] = -s[i]
    return s

def tmpr(num):
    return 10**10*0.8**np.floor(num/300)

def prepartGen(input_seq):
    n = len(input_seq)
    p = [random.randint(0, n-1) for _ in range(n)]
    return p

def p_trans_res(input_seq, p):
    n = len(input_seq)
    a_ = [0 for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if p[j] == i:
                a_[i] += input_seq[j]
    return KK(a_)

def random_neighbors_p(partition):
    p = copy.deepcopy(partition)
    n = len(partition)
    for t in range(n):
        i, j = random.sample(range(n), 2)
        if p[i] != j:
            p[i] = j
            break
    return p


# In[67]:


def repeated_random(input_seq, num_iter):
    n = len(input_seq)
    s = random_solutionGen(n)
    for i in range(num_iter):
        s_ = random_solutionGen(n)
        if residue(input_seq, s_) < residue(input_seq, s):
            s = s_
    return residue(input_seq, s)

def hill_climb(input_seq, num_iter):
    n = len(input_seq)
    s = random_solutionGen(n)
    for i in range(num_iter):
        s_ = random_neighbors(s)
        if residue(input_seq, s_) < residue(input_seq, s):
            s = s_
    return residue(input_seq, s)

def simulated_anneal(input_seq, num_iter):
    n = len(input_seq)
    s = random_solutionGen(n)
    s__ = copy.deepcopy(s)
    for i in range(num_iter):
        s_ = random_neighbors(s)
        if residue(input_seq, s_) < residue(input_seq, s):
            s = s_
        else:
            if random.random() < np.exp(-(residue(input_seq, s_)-residue(input_seq, s))/tmpr(i)):
                s = s_
        if residue(input_seq, s) < residue(input_seq, s__):
            s__ = s
    return residue(input_seq, s__)

def prepart_rand(input_seq, num_iter):
    partition = prepartGen(input_seq)
    for i in range(num_iter):
        partition_ = prepartGen(input_seq)
        if p_trans_res(input_seq, partition_) < p_trans_res(input_seq, partition):
            partition = partition_
    return p_trans_res(input_seq, partition)

def prepart_hill(input_seq, num_iter):
    p = prepartGen(input_seq)
    for i in range(num_iter):
        p_ = random_neighbors_p(p)
        if p_trans_res(input_seq, p_) < p_trans_res(input_seq, p):
            p = p_
    return p_trans_res(input_seq, p)

def prepart_simul(input_seq, num_iter):
    p = prepartGen(input_seq)
    p__ = copy.deepcopy(p)
    for i in range(num_iter):
        p_ = random_neighbors_p(p)
        if p_trans_res(input_seq, p_) < p_trans_res(input_seq, p):
            p = p_
        else:
            if random.random() < np.exp(-(p_trans_res(input_seq, p_)-p_trans_res(input_seq, p))/tmpr(i)):
                p = p_
        if p_trans_res(input_seq, p) < p_trans_res(input_seq, p__):
            p__ = p
    return p_trans_res(input_seq, p__)


# In[46]:


def main_fc(flag, algr, input_file, num_iter):
    indata = read_file(input_file)
    
    if flag !=0:
        import time
        start = time.time()
        
    if algr == 0:
        res = KK(indata)
    if algr == 1:
        res = repeated_random(indata, num_iter)
    if algr == 2:
        res = hill_climb(indata, num_iter)
    if algr == 3:
        res = simulated_anneal(indata, num_iter)
    if algr == 11:
        res = prepart_rand(indata, num_iter)
    if algr == 12:
        res = prepart_hill(indata, num_iter)
    if algr == 13:
        res = prepart_simul(indata, num_iter)
        
    if flag !=0:
        end = time.time()
        print('The running time per trial is: {0} ms'.format((end-start) * 1000/numtrials))
    return res


# In[36]:


if __name__ == "__main__":
    
    import sys
    import heapq
    import numpy as np
    import random
    import copy
    
    num_iter = 3000
    output = main_fc(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], num_iter)
    print(output)

