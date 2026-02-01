import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from math import *
import copy

OUTPUT_EPSILON = '{:0.8f}'
ZEIDEL_EPSILON = 1e-6
WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def find_leading(matrix, b, step) -> tuple[list[list[float]], list[float]]:
        leading_row = step
        for i in range(step+1, len(matrix)):
            if abs(matrix[i][step]) > abs(matrix[leading_row][step]):
                leading_row = i
        matrix[leading_row], matrix[step] = matrix[step], matrix[leading_row]
        b[leading_row], b[step] = b[step], b[leading_row]
        return matrix, b


def gauss(matrix, b) -> tuple[list[list[float]], list[float], list[float]]:
    for step in range(len(matrix)-1):
        matrix, b = find_leading(copy.deepcopy(matrix), b, step)
        for i in range(len(matrix)-1, step, -1):
            division = matrix[i][step] / matrix[i-1][step]
            for j in range(len(matrix)-1, step-1, -1): 
                matrix[i][j] -= matrix[i-1][j] * division
            b[i] -= b[i-1] * division
            matrix[i][step] = 0

    x = [0] * len(matrix)
    for step in range(len(matrix)-1, -1, -1):
        x[step] = (b[step] - sum([matrix[step][i] * x[i] for i in range(len(matrix))])) / matrix[step][step]
    return matrix, b, x



def zeidel(matrix, b) -> tuple[bool, list[list[float]], list[float], list[float], int]:
    for step in range(len(matrix)-1):
        matrix, b = find_leading(copy.deepcopy(matrix), b, step)
    print(matrix)
    checked = all(abs(matrix[i][i]) > sum(abs(matrix[i][j]) for j in range(len(matrix)) if i != j) for i in range(len(matrix)))
    print(f"Usage of Zeidel is possible:\t{checked}")
    if not checked:
        return checked, matrix, b, None, None
    counter = 0
    x = [0] * len(matrix)
    x_prev = [inf] * len(matrix)
    while all([abs(x[i] - x_prev[i]) > ZEIDEL_EPSILON for i in range(len(matrix))]):
    # for _ in range(10):
        x_prev = copy.deepcopy(x)
        for i in range(len(matrix)):
            print(x[i], b[i], matrix[i][i], sum([matrix[i][j] * x[j] for j in range(len(matrix)) if i != j]))
            x[i] = (b[i] - sum([matrix[i][j] * x[j] for j in range(len(matrix)) if i != j])) / matrix[i][i]
        print(x)
        counter += 1
    return checked, matrix, b, x, counter



if __name__ == '__main__':
    A = [
            [  6.1,    -2.2,    -1.2,    -3.3  ],
            [  7.2,     0.9,     1.8,    -4.1  ],
            [  2.8,     3.3,     1.1,     2.5  ],
            [ -1.5,     1.0,     6.3,     0.8  ]
        ]
    b     = [ -0.4,     3.1,    10.8,    -5.0  ]
    x_true= [  1.0,     2.0,    -1.0,     1.0  ]

    A_gauss, b_gauss, ans_gauss = gauss(copy.deepcopy(A), copy.deepcopy(b))
    with open(fr'{WORKDIR}\Results\Guass.csv', 'w') as file:
        file.writelines('\n'.join(['\t'.join([OUTPUT_EPSILON.format(item) for item in A_gauss[i]]) + '\t|\t' + OUTPUT_EPSILON.format(b_gauss[i]) for i in range(len(A))]))
        file.write('\n' * 3)
        file.write('\n'.join([OUTPUT_EPSILON.format(item) for item in ans_gauss]))
        file.write('\n' * 3)
        file.write(str(max(abs(np.array(ans_gauss) - np.array(x_true)))))
    
    zeidel_possible, A_zeidel, b_zeidel, ans_zeidel, num_iter = zeidel(copy.deepcopy(A), copy.deepcopy(b))
    with open(fr'{WORKDIR}\Results\Zeidel.csv', 'w') as file:
        file.write('Usage of Zeidel is possible:' + '\t' + str(zeidel_possible) + '\n' * 3)
        file.writelines('\n'.join(['\t'.join([OUTPUT_EPSILON.format(item) for item in A_zeidel[i]]) + '\t|\t' + OUTPUT_EPSILON.format(b_zeidel[i]) for i in range(len(A))]))
        if not zeidel_possible:
            exit(0)
        file.write('\n' * 3)
        file.write('\n'.join([OUTPUT_EPSILON.format(item) for item in ans_zeidel]))
        file.write('\n' * 3)
        file.write(str(num_iter))
        file.write('\n' * 3)
        file.write(str(max(abs(np.array(ans_zeidel) - np.array(x_true)))))