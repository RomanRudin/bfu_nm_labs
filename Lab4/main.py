import numpy as np
from pathlib import Path
from math import *
from typing import Callable
import pandas as pd

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def interpolate(y: Callable, yn: Callable, x: float, xk: list[float]) -> pd.DataFrame:
    L = lambda x: sum(y[i] * prod([(x - xi) for xi in xk]) / prod([xk[i] - xj for xj in xk if xj != x[i]]) for i in range(len(xk)))
    max_yn_x = np.linspace(xk[0], xk[-1], int(abs(xk[-1] - xk[0]) * 1000))
    M = max(abs(yn(max_yn_x)))
    R_max = lambda x: M / factorial(len(xk) + 1) * abs(prod([x - xi for xi in xk]))
    table = {'xk': 'yk'}
    table.update({OUTPUT_EPSILON(xk[i]): OUTPUT_EPSILON(y(xk[i])) for i in range(len(xk))})
    df = pd.DataFrame(table, index=[0])
    with open(fr'{WORKDIR}\Results\Interpolation.csv', 'w') as file:
        df.to_csv(path_or_buf=file, sep='\t', lineterminator='\n', index=False)
        file.write('\n' * 2 + 'Method error' + '\n')
        

    #TODO

def differentiate(f: Callable, n: int, a: float, b: float, f1: Callable, f2: Callable) -> pd.DataFrame:
    def calculate_error(x: list[float], x_actual: list[float]) -> float:
        return max([abs(xi - xi_actual)] for xi, xi_actual in zip(x, x_actual))
    
    h = (b - a) / (n + 1)
    x = [a + i * h for i in range(1, n+1)]
    y = [f(i) for i in [a] + x + [b]]
    f_left, f_right, f_central, f_second, f1_actual, f2_actual = [], [], [], [], [], [] 
    for i in range (1, n + 1):
        f_left.append((y[i] - y[i-1]) / h)
        f_right.append((y[i+1] - y[i]) / h)
        f_central.append((y[i+1] - y[i-1]) / (2 * h))
        f_second.append((y[i-1] - 2 * y[i] + y[i+1]) / h**2)
        f1_actual.append(f1(x[i-1]))
        f2_actual.append(f2(x[i-1]))
    errors_df = pd.DataFrame({
            'Left': calculate_error(f_left, f1_actual),
            'Right': calculate_error(f_right, f1_actual),
            'Central': calculate_error(f_central, f1_actual),
            'Second': calculate_error(f_second, f2_actual),
        })
    answer_df = pd.DataFrame({
            'x': x,
            'Left': f_left,
            'Right': f_right,
            'Central': f_central,
            'Second': f_second,
            'Actual f\'(x)': f1_actual,
            'Actual f\'\'(x)': f2_actual,
        })
    with open(fr'{WORKDIR}\Results\Differentiation.csv', 'w') as file:
        # answer_df.to_csv(path_or_buf=file, sep='\t', lineterminator='\n', float_format=OUTPUT_EPSILON_DF)
        names = list(answer_df.columns)
        max_sym = max([len(i) for i in names + [OUTPUT_EPSILON(-10)]])
        file.write(' | '.join([name.rjust(max_sym) for name in names]) + '\n')
        for _, row in answer_df.iterrows():
            file.write(' | '.join([OUTPUT_EPSILON(item).rjust(max_sym) for item in row]) + '\n')
        file.write('\n' * 2 + 'Error of each method:' + '\n')
        errors_df.to_csv(path_or_buf=file, sep='\t', lineterminator='\n', float_format=OUTPUT_EPSILON_DF, index=False)
    #TODO Make error calculation?
    return answer_df, errors_df

    
if __name__ == "__main__":
    OUTPUT_EPSILON = '{:0.4f}'.format
    y = lambda x: cos(2*x)
    y5 = lambda x: -32 * np.sin(2*x)
    x = 1.18
    xk = [1.00, 1.10, 1.20, 1.30]
    interpolate(y, y5, x, xk)

    OUTPUT_EPSILON = '{:0.4f}'.format
    OUTPUT_EPSILON_DF = '%.4f'
    n = 5
    f = lambda x: cos(2*x) ** 2
    f1 = lambda x: -4 * cos(2*x) *  sin(2*x)
    f2 = lambda x: 8 * (sin(2*x) ** 2 - cos(2*x) ** 2)
    a, b = 0, 1
    differentiate(f, n, a, b, f1, f2)