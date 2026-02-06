import numpy as np
from pathlib import Path
from math import *
from typing import Callable

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def logger(func_name: str, func: Callable) -> None:
    func_name = func_name.lower()
    I_prox, proximity, delta, n = func(f, INTERVAL[0], INTERVAL[1], f2, f4)
    with open(fr'{WORKDIR}\Results\{func_name.title().replace(" ", "")}.csv', 'w') as file:
        file.write(f'{func_name.upper()} METHOD' + '\n' * 2)
        file.write(f'Exact value: \t {OUTPUT_EPSILON(I)}' + '\n')
        file.write(f'Obtained value: \t {OUTPUT_EPSILON(I_prox)}' + '\n')
        file.write(f'Calculation error: \t {OUTPUT_EPSILON(proximity)}' + '\n')
        file.write(f'Relative error: \t {OUTPUT_EPSILON((abs(I-I_prox)/I) * 100)}%' + '\n')
        file.write(f'Method error: \t {OUTPUT_EPSILON(delta)}' + '\n' * 2)
        file.write(f'Number of iterations: \t {int(log2(n))}' + '\n')
        file.write(f'n: \t {(n)}' + '\n')



def left_rectangles(f: Callable, a: float, b: float, f2: Callable, f4: Callable) -> float:
    n = 2
    I_prev = -inf
    I = inf
    while abs(I - I_prev) > EPSILON:
        n *= 2
        I_prev = I
        h = (b - a) / n
        I = h * sum([f(a + i * h) for i in range(0, n)])
    proximity = I - I_prev
    x = np.linspace(a, b, (b-a)*1000)
    delta = max(abs(f2(x))) / 24 * (b-a) * h**2
    return I, proximity, delta, n

def right_rectangles(f: Callable, a: float, b: float, f2: Callable, f4: Callable) -> float:
    n = 2
    I_prev = -inf
    I = inf
    while abs(I - I_prev) > EPSILON:
        n *= 2
        I_prev = I
        h = (b - a) / n
        I = h * sum([f(a + i * h) for i in range(1, n+1)])
    proximity = I - I_prev
    x = np.linspace(a, b, (b-a)*1000)
    delta = max(abs(f2(x))) / 24 * (b-a) * h**2
    return I, proximity, delta, n

def middle_rectangles(f: Callable, a: float, b: float, f2: Callable, f4: Callable) -> float:
    n = 2
    I_prev = -inf
    I = inf
    while abs(I - I_prev) > EPSILON:
        n *= 2
        I_prev = I
        h = (b - a) / n
        I = h * sum([f(a + h/2 + i * h) for i in range(0, n)])
    proximity = I - I_prev
    x = np.linspace(a, b, (b-a)*1000)
    delta = max(abs(f2(x))) / 24 * (b-a) * h**2
    return I, proximity, delta, n



def trapezoids(f: Callable, a: float, b: float, f2: Callable, f4: Callable) -> float:
    n = 2
    I_prev = -inf
    I = inf
    while abs(I - I_prev) > EPSILON:
        n *= 2
        I_prev = I
        h = (b - a) / n
        I = h * (f(a) + f(b) / 2 + sum([f(a + i * h) for i in range(1, n-1)]))
    proximity = I - I_prev
    x = np.linspace(a, b, (b-a)*1000)
    delta = max(abs(f2(x))) / 12 * (b-a) * h**2
    return I, proximity, delta, n




def simpson(f: Callable, a: float, b: float, f2: Callable, f4: Callable) -> float:
    n = 2
    I_prev = -inf
    I = inf
    while abs(I - I_prev) > EPSILON:
        n *= 2
        I_prev = I
        h = (b - a) / n
        I = h / 3 * (f(a) + f(b) + 4 * sum([f(a + i * h) for i in range(1, n, 2)]) + 2 * sum([f(a + i * h) for i in range(2, n, 2)]))
    proximity = I - I_prev
    x = np.linspace(a, b, (b-a)*1000)
    delta = (b-a) / 2880 * h**4 * max(abs(f4(x)))
    return I, proximity, delta, n




if __name__ == "__main__":
    f = lambda x: (1 + x)**(1/2)
    EPSILON = 1e-4
    OUTPUT_EPSILON = '{:0.6f}'.format
    INTERVAL = [0, 1]
    F = lambda x: (2/3)*((1 + x)**3)**(1/2)
    f2 = lambda x: - 1 / (4 * (1 + x)**(3/2))
    f4 = lambda x: - 15 / (16 * (1 + x)**(7/2))

    I = F(INTERVAL[1]) - F(INTERVAL[0])
    logger('left rectangles', left_rectangles)
    logger('right rectangles', right_rectangles)
    logger('middle rectangles', middle_rectangles)
    logger('trapezoids', trapezoids)
    logger('simpson', simpson)