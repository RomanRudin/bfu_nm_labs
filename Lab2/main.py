from pathlib import Path
from math import *
import copy

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def calculate_row_condition(row: list[float], step: int) -> float:
    return abs(row[step]) - sum(abs(row[i]) for i in range(len(row)) if i != step)


def change_to_leading(matrix: list[list[float]], b: list[float], step: int) -> None:
        leading_row = step
        for i in range(step+1, len(matrix)):
            if abs(matrix[i][step]) > abs(matrix[leading_row][step]):
                leading_row = i
        matrix[leading_row], matrix[step] = matrix[step], matrix[leading_row]
        b[leading_row], b[step] = b[step], b[leading_row]



def gauss(matrix: list[list[float]], b: list[float]) -> tuple[list[list[float]], list[float], list[float]]:
    for step in range(len(matrix)-1):
        change_to_leading(matrix, b, step)
        for i in range(len(matrix)-1, step, -1):
            division = matrix[i][step] / matrix[i-1][step]
            for j in range(len(matrix)-1, step-1, -1): 
                matrix[i][j] -= matrix[i-1][j] * division
            b[i] -= b[i-1] * division
            matrix[i][step] = 0

    x = [0] * len(matrix)
    for step in range(len(matrix)-1, -1, -1):
        x[step] = (b[step] - sum([matrix[step][i] * x[i] for i in range(len(matrix))])) / matrix[step][step]
    return matrix, b, x, 0



def zeidel(matrix: list[list[float]], b: list[float]) -> tuple[bool, list[list[float]], list[float], list[float], int]:
    for step in range(len(matrix)-1):
        change_to_leading(matrix, b, step)
    checked = all(calculate_row_condition(matrix[i], i) > 0 for i in range(len(matrix)))
    counter = 0
    x = [0] * len(matrix)
    x_prev = [inf] * len(matrix)
    error = inf
    while error > ZEIDEL_EPSILON:
    # for _ in range(3):
        x_prev = copy.deepcopy(x)
        for i in range(len(matrix)):
            x[i] = (b[i] - sum([matrix[i][j] * x[j] for j in range(len(matrix)) if i != j])) / matrix[i][i]
        counter += 1
        error = max([abs(x[i] - x_prev[i]) for i in range(len(matrix))]) 
        print([OUTPUT_EPSILON.format(item) for item in x], OUTPUT_EPSILON.format(error))
    return checked, matrix, b, x, counter, error



if __name__ == '__main__':
    A = [
            [  0.5,     8.0,    -1.2,    -1.3  ],
            [  7.2,     0.9,     1.8,    -2.2  ],
            [  0.8,     1.1,     1.1,     5.3  ],
            [ -1.5,     1.0,     6.3,     0.8  ]
        ]
    b     = [ 16.4,     5.0,     7.2,    -5.0  ]
    x_true= [  1.0,     2.0,    -1.0,     1.0  ]

    eps = len(str(A[0][0]).split('.')[-1])
    print(f'Please, change OUTPUT_EPSILON to {eps + 1} or {eps + 2} and ZEIDEL_EPSILON {eps}. Change ZEIDEL_ERROR_OUTPUT_EPSILON if you need too.')
    OUTPUT_EPSILON = '{:0.2f}'
    ZEIDEL_EPSILON = 1e-1
    ZEIDEL_ERROR_OUTPUT_EPSILON = '{:0.3f}'

    A_gauss, b_gauss, ans_gauss, error = gauss(copy.deepcopy(A), copy.deepcopy(b))
    with open(fr'{WORKDIR}\Results\Guass.csv', 'w') as file:
        file.write('Gauss method:' + '\n' * 2)
        file.write('Prepared matrix:\n')
        file.writelines('\n'.join(['\t'.join([OUTPUT_EPSILON.format(item) for item in A_gauss[i]]) + '\t|\t' + OUTPUT_EPSILON.format(b_gauss[i]) for i in range(len(A))]))
        file.write('\n' * 3)
        file.write('Vector of approximate solution:\n')
        file.write('\n'.join([OUTPUT_EPSILON.format(item) for item in ans_gauss]))
        file.write('\n' * 3)
        file.write('Calculation error:\n')
        file.write(f'{error}. The method is an exact method')
    
    zeidel_possible, A_zeidel, b_zeidel, ans_zeidel, num_iter, error = zeidel(copy.deepcopy(A), copy.deepcopy(b))
    with open(fr'{WORKDIR}\Results\Zeidel.csv', 'w') as file:
        file.write('Iterative Gauss-Zeidel method:' + '\n' * 2)
        file.write('Usage of Zeidel is possible:' + '\t' + str(zeidel_possible) + '\n' * 3)
        file.write('Prepared matrix:\n')
        file.writelines('\n'.join(['\t'.join([OUTPUT_EPSILON.format(item) for item in A_zeidel[i]]) + '\t|\t' + OUTPUT_EPSILON.format(b_zeidel[i]) for i in range(len(A))]))
        if not zeidel_possible:
            file.write('\n' * 3)
            file.write('A sufficient (but not necessary!) condition is not met. Perhaps a matrix transformation is necessary.\n')
            file.write('Look at page 301 at Demidovich-Maron')
        file.write('\n' * 3)
        file.write('Vector of approximate solution:\n')
        file.write('\n'.join([OUTPUT_EPSILON.format(item) for item in ans_zeidel]))
        file.write('\n' * 3)
        file.write('Number of iterations:\n')
        file.write(str(num_iter))
        file.write('\n' * 3)
        file.write('Calculated error:\n')
        file.write(ZEIDEL_ERROR_OUTPUT_EPSILON.format(error))