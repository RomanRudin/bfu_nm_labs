import numpy as np
from pathlib import Path
from math import *
from typing import Callable
import pandas as pd

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def interpolate(y: Callable, yn: Callable, x: float, xk: list[float]) -> None:
    L = lambda x: sum(y[i] * prod([(x - xi) for xi in xk]) / prod([xk[i] - xj for xj in xk if xj != x[i]]) for i in range(len(xk)))
    max_yn_x = np.linspace(xk[0], xk[-1], int(abs(xk[-1] - xk[0]) * 1000))
    M = max(abs(yn(max_yn_x)))
    R_max = lambda x: M / factorial(len(xk) + 1) * abs(prod([x - xi for xi in xk]))
    df = pd.DataFrame(
        {
            'x': xk,
            'y': [y(i) for i in xk],
        }
    )
    #TODO

def differentiate(f: Callable, n: int, a: float, b: float, f1: Callable, f2: Callable) -> pd.DataFrame:
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
    df = pd.DataFrame(
        {
            'x': x,
            'Left': f_left,
            'Right': f_right,
            'Central': f_central,
            'Second': f_second,
            'Actual f\'(x)': f1_actual,
            'Actual f\'\'(x)': f2_actual,
        }
    )
    with open(fr'{WORKDIR}\Results\Differentiation.csv', 'w') as file:
        df.to_csv(path_or_buf=file, sep='\t', lineterminator='\n', float_format=OUTPUT_EPSILON_DF)
        #TODO Make error calculation
    return df

    
if __name__ == "__main__":
    OUTPUT_EPSILON = '{:0.2f}'
    y = lambda x: cos(2*x)
    y5 = lambda x: -32 * np.sin(2*x)
    x = 1.18
    xk = [1.00, 1.10, 1.20, 1.30]
    interpolate(y, y5, x, xk)

    OUTPUT_EPSILON = '{:0.5f}'
    OUTPUT_EPSILON_DF = '%.5f'
    n = 5
    f = lambda x: cos(2*x) ** 2
    f1 = lambda x: -4 * cos(2*x) *  sin(2*x)
    f2 = lambda x: 8 * (sin(2*x) ** 2 - cos(2*x) ** 2)
    a, b = 0, 1
    differentiate(f, n, a, b, f1, f2)