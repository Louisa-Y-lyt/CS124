#!/usr/bin/env python
# coding: utf-8

# In[23]:


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math

def strassen(a, b, n0):
    """
    Multiply two matrices a and b together using Strassen's algorithm
    """
    # Get dimensions of input matrices
    a_rows, a_cols = len(a), len(a[0])
    b_rows, b_cols = len(b), len(b[0])
    
    # Check that the matrices are compatible for multiplication
    if a_cols != b_rows:
        raise ValueError("Matrices are not compatible for multiplication")
    
    # Pad matrices with zeros to make their dimensions powers of two
    max_dim = max(a_rows, a_cols, b_cols)
    padded_dim = 2 ** math.ceil(math.log2(max_dim))
    a_padded = [[0 for _ in range(padded_dim)] for _ in range(padded_dim)]
    b_padded = [[0 for _ in range(padded_dim)] for _ in range(padded_dim)]
    for i in range(a_rows):
        for j in range(a_cols):
            a_padded[i][j] = a[i][j]
    for i in range(b_rows):
        for j in range(b_cols):
            b_padded[i][j] = b[i][j]

    # Perform Strassen's algorithm
    result_padded = _strassen_helper(a_padded, b_padded, n0)

    # Remove padding from result matrix
    result = [[0 for _ in range(b_cols)] for _ in range(a_rows)]
    for i in range(a_rows):
        for j in range(b_cols):
            result[i][j] = result_padded[i][j]
    
    return result

def _strassen_helper(a, b, n0):
    """
    Helper function to perform Strassen's algorithm on padded matrices
    """
    # Get dimensions of input matrices
    dim = len(a)
    
    # Base case: if dimensions are 1x1, return their product
    if dim <= n0:
        return simple_mulMat(a, b)
    
    # Split matrices into quadrants
    half_dim = dim // 2
    a11 = [a[i][:half_dim] for i in range(half_dim)]
    a12 = [a[i][half_dim:] for i in range(half_dim)]
    a21 = [a[i][:half_dim] for i in range(half_dim, dim)]
    a22 = [a[i][half_dim:] for i in range(half_dim, dim)]
    b11 = [b[i][:half_dim] for i in range(half_dim)]
    b12 = [b[i][half_dim:] for i in range(half_dim)]
    b21 = [b[i][:half_dim] for i in range(half_dim, dim)]
    b22 = [b[i][half_dim:] for i in range(half_dim, dim)]
    
    # Compute products of submatrices
    p1 = _strassen_helper(a11, _sub_matrices(b12, b22),n0)
    p2 = _strassen_helper(_add_matrices(a11, a12), b22,n0)
    p3 = _strassen_helper(_add_matrices(a21, a22), b11,n0)
    p4 = _strassen_helper(a22, _sub_matrices(b21, b11), n0)
    p5 = _strassen_helper(_add_matrices(a11, a22), _add_matrices(b11, b22),n0)
    p6 = _strassen_helper(_sub_matrices(a12, a22), _add_matrices(b21, b22),n0)
    p7 = _strassen_helper(_sub_matrices(a21, a11), _add_matrices(b11, b12),n0)

    r11 = _add_matrices(_add_matrices(_sub_matrices(p4, p2),p5),p6)
    r12 = _add_matrices(p1,p2)
    r21 = _add_matrices(p3,p4)
    r22 = _add_matrices(_add_matrices(_sub_matrices(p1, p3),p5),p7)
    
    # Grouping the results obtained in a single matrix:
    R = [[0 for j in range(0, dim)] for i in range(0, dim)]
    for i in range(0, half_dim):
        for j in range(0, half_dim):
            R[i][j] = r11[i][j]
            R[i][j + half_dim] = r12[i][j]
            R[i + half_dim][j] = r21[i][j]
            R[i + half_dim][j + half_dim] = r22[i][j]
    return R
    
def _sub_matrices(a,b):
    n = len(a)
    rslt = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            rslt[i][j] = a[i][j]-b[i][j]
    return rslt

def _add_matrices(a,b):
    n = len(a)
    rslt = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            rslt[i][j] = a[i][j]+b[i][j]
    return rslt

def simple_mulMat(mat1, mat2):
    # List to store matrix multiplication result
    n = len(mat1)
    rslt = [[0 for _ in range(n)] for _ in range(n)]
  
    for i in range(n):
        for j in range(n):
            for k in range(n):
                rslt[i][j] += mat1[i][k] * mat2[k][j]
    return rslt


# In[ ]:


def read_file(input_file, n):
    with open(input_file, 'r') as f:
        rawdata = np.genfromtxt(f, dtype=int, delimiter=',')
    mat1 = [[0 for _ in range(n)] for _ in range(n)]
    mat2 = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        mat1[i][:n] = rawdata[i*n:(i+1)*n]
    for i in range(n,2*n):
        mat2[i-n][:n] = rawdata[i*n:(i+1)*n]
    return mat1, mat2

def main_fc(flag, d, input_file):
    mat1, mat2 = read_file(input_file, d)
    if flag !=0:
        import time
        start = time.time()
    mat3 = strassen(mat1, mat2, 32)
    res = np.diagonal(mat3)
    if flag !=0:
        end = time.time()
        print('The running time per trial is: {0} ms'.format((end-start) * 1000/numtrials))
    return res


# In[8]:


if __name__ == "__main__":
    
    import sys
    import math
    import numpy as np
    
    output = main_fc(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
    print('\n'.join(str(v) for v in output), )


# In[ ]:




