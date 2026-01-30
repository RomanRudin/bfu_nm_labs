import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from math import *

def Gauss(A, b):
    def find_leading(A, step):
        leading_row = step
        for i in range(step+1, len(A)):
            if abs(A[i][step]) > abs(A[leading_row][step]):
                leading_row = i
        A[leading_row], A[step] = A[step], A[leading_row]
        b[leading_row], b[step] = b[step], b[leading_row]
        # print(leading_row)
        # print(b)
        return A, b
    
    for step in range(len(A)-1):
        A, b = find_leading(A, step)
        for i in range(len(A)-1, step, -1):
            division = A[i][step] / A[i-1][step]
            for j in range(len(A)-1, step-1, -1): 
                # print(A[i][j],division, A[i-1][j], sep='\t')
                A[i][j] -= A[i-1][j] * division
            b[i] -= b[i-1] * division
            A[i][step] = 0
            # print(b)
            # print('\n'.join([''.join(['{:8.4f}'.format(item) for item in row]) + '\t' + str(b[k]) for k, row in enumerate(A)]), end='\n'*2) # print(A, end='\n'*2)

    ans = [0] * len(A)
    for step in range(len(A)-1, -1, -1):
        # print(step, ans, b[step])
        ans[step] = (b[step] - sum([A[step][i] * ans[i] for i in range(len(A))])) / A[step][step]
        # print(ans[step])
    return A, ans


def Zeidel(A, b):
    pass    



if __name__ == '__main__':
    A = [
            [  6.1,    -2.2,    -1.2,    -3.3  ],
            [  7.2,     0.9,     1.8,    -4.1  ],
            [  2.8,     3.3,     1.1,     2.5  ],
            [ -1.5,     1.0,     6.3,     0.8  ]
        ]
    b     = [ -0.4,     3.1,    10.8,    -5.0  ]
    x_true= [  1.0,     2.0,    -1.0,     1.0  ]

    A_gauss, ans_gauss = Gauss(A, b) 
    print('\n'.join([''.join(['{:8.4f}'.format(item) for item in row]) for row in A_gauss]))
    print('\t'.join([str(item) for item in ans_gauss]))
    print(Zeidel(A, b))