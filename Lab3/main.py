import numpy as np
import pandas as pd
from pathlib import Path
from math import *

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))



def simpson(f: lambda x: float, a: float, b: float, n: int) -> float:
    h = (b - a) / n
    I = h / 3 * (f(a) + f(b) + 4 * sum([f(a + i * h) for i in range(1, n, 2)]) + 2 * sum([f(a + i * h) for i in range(2, n, 2)]))
    delta = (b - a) / 2880 * h**4 * max(abs())



def left_rectangles(f: lambda x: float, a: float, b: float, n: int, f2: lambda x: float) -> float:
    h = (b - a) / n
    I = h * sum([f(a + i * h) for i in range(0, n)])
    delta = max(abs(f2))*(b-a)*h**2
    return I, delta

def right_rectangles(f: lambda x: float, a: float, b: float, n: int, f2: lambda x: float) -> float:
    h = (b - a) / n
    I = h * sum([f(a + i * h) for i in range(1, n+1)])
    delta = max(abs(f2))*(b-a)*h**2
    return I, delta

def medium_rectangles(f: lambda x: float, a: float, b: float, n: int, f2: lambda x: float) -> float:
    h = (b - a) / n
    I = h * sum([f(a + h/2 + i * h) for i in range(0, n)])
    delta = max(abs(f2))*(b-a)*h**2
    return I, delta



def trapezoid(f: lambda x: float, a: float, b: float, n: int) -> float:
    h = (b - a) / n
    I = h * (f(a) + f(b)) / 2 + sum([f(a + i * h) for i in range(0, n)])
    delta = max(abs(f2))*(b-a)*h**2
    pass



if __name__ == "__main__":
    f = lambda x: (1 + x)**(1/2)
    EPSILON = 1e-4
    INTERVAL = [0, 1]
    F = lambda x: (2/3)*((1 + x)**3)**(1/2)
    f2 = 0
    f4 = 0