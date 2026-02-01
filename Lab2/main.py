import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from math import *

OUTPUT_EPSILON = '{:0.20f}'
WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def Gauss(A, b):
    def find_leading(A, step):
        leading_row = step
        for i in range(step+1, len(A)):
            if abs(A[i][step]) > abs(A[leading_row][step]):
                leading_row = i
        A[leading_row], A[step] = A[step], A[leading_row]
        b[leading_row], b[step] = b[step], b[leading_row]
        return A, b
    
    for step in range(len(A)-1):
        A, b = find_leading(A, step)
        for i in range(len(A)-1, step, -1):
            division = A[i][step] / A[i-1][step]
            for j in range(len(A)-1, step-1, -1): 
                A[i][j] -= A[i-1][j] * division
            b[i] -= b[i-1] * division
            A[i][step] = 0

    ans = [0] * len(A)
    for step in range(len(A)-1, -1, -1):
        ans[step] = (b[step] - sum([A[step][i] * ans[i] for i in range(len(A))])) / A[step][step]
    return A, ans



def Zeidel(A, b) -> list[int]:
    functions_array = []
    for i in range(len(A)):
        functions_array.append(lambda x: (b[i] - sum([A[i][j] * x[j] for j in range(len(A)) if i != j]) ) / A[i][i])
    x = [0] * len(A)
    x_prev = [inf] * len(A)
    while all([abs(x[i] - x_prev[i]) > 1e-6 for i in range(len(A))]):
    # for _ in range(4):
        for i in range(len(A)): x_prev[i] = x[i]
        for i in range(len(A)):
            x[i] = functions_array[i](x)
    return x




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
    with open(fr'{WORKDIR}\Results\Guass.csv', 'w') as file:
        file.write('\n'.join(['\t'.join([OUTPUT_EPSILON.format(item) for item in row]) for row in A_gauss]))
        file.write('\n' * 3)
        file.write('\n'.join([OUTPUT_EPSILON.format(item) for item in ans_gauss]))
    
    ans_zeidel = Zeidel(A, b)
    with open(fr'{WORKDIR}\Results\Zeidel.csv', 'w') as file:
        file.write('\n'.join([OUTPUT_EPSILON.format(item) for item in ans_zeidel]))