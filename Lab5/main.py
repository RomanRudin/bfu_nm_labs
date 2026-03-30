import numpy as np
from pathlib import Path
import pandas as pd
from typing import Callable
from math import log10

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def f1(x, y) -> float:
    return np.cos(1.5 * x + y) + (x - y)


def f2_system(x, Y) -> np.ndarray:
    y1, y2 = Y
    dy1 = y2
    dy2 = np.cos(1.5 + x) + 0.1 * y1**2
    return np.array([dy1, dy2])


def euler_cauchy(f: Callable, n: int) -> tuple:
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1)

    for i in range(n):
        k1 = f(x[i], y[i])
        y_pred = y[i] + h * k1
        k2 = f(x[i] + h, y_pred)
        y[i + 1] = y[i] + h * (k1 + k2) / 2

    return x, y


def runge_kutta4(f: Callable, n: int) -> tuple:
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1)

    for i in range(n):
        k1 = h * f(x[i], y[i])
        k2 = h * f(x[i] + h/2, y[i] + k1/2)
        k3 = h * f(x[i] + h/2, y[i] + k2/2)
        k4 = h * f(x[i] + h, y[i] + k3)

        y[i + 1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6

    return x, y


def rk4_system(f: Callable, n: int) -> tuple:
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    Y = np.zeros((n + 1, 2))
    Y[0] = [0, 1]

    for i in range(n):
        k1 = h * f(x[i], Y[i])
        k2 = h * f(x[i] + h/2, Y[i] + k1/2)
        k3 = h * f(x[i] + h/2, Y[i] + k2/2)
        k4 = h * f(x[i] + h, Y[i] + k3)

        Y[i + 1] = Y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6

    return x, Y


def adams3(f: Callable, n: int) -> tuple:
    h = (b - a) / n
    x, Y = rk4_system(f, n)

    for i in range(2, n):
        f_i  = f(x[i], Y[i])
        f_i1 = f(x[i-1], Y[i-1])
        f_i2 = f(x[i-2], Y[i-2])

        Y[i+1] = Y[i] + h * (23*f_i - 16*f_i1 + 5*f_i2) / 12

    return x, Y[:, 0]


def adams4(f: Callable, n: int) -> tuple:
    h = (b - a) / n
    x, Y = rk4_system(f, n)

    for i in range(3, n):
        f_i  = f(x[i], Y[i])
        f_i1 = f(x[i-1], Y[i-1])
        f_i2 = f(x[i-2], Y[i-2])
        f_i3 = f(x[i-3], Y[i-3])

        Y[i+1] = Y[i] + h * (55*f_i - 59*f_i1 + 37*f_i2 - 9*f_i3) / 24

    return x, Y[:, 0]


def calculate(method: Callable, f: Callable) -> tuple:
    n = 2

    while True:
        x_prev, y_prev = method(f, n-1)
        x_last, y_last = method(f, 2 * n-1)
        y_last_match = y_last[::2]
        diff = np.max(np.abs(y_prev - y_last_match))
        if diff < EPSILON:
            print(diff)
            print(f"n:{2*n}")
            return n, x_prev, y_prev, 2*n, x_last, y_last
        n *= 2


def save(method_name: str, n_prev: int, x_prev: list, y_prev: list, n_last: int, x_last: list, y_last: list) -> None:
    y_prev_sel = list(map(lambda x: str(x), y_prev[-8:]))
    for i in range(8):
        y_prev_sel.insert(2*i, '-')
    y_last_sel = y_last[-16:]

    diff = [abs(last - float(prev)) if prev != '-' else '-' for last, prev in zip(y_last_sel, y_prev_sel)]

    df = pd.DataFrame({
        "x_k": x_last[-16:],
        "y_prev": y_prev_sel,
        "y_last": y_last_sel,
        "difference": diff
    })

    df = df.round(int(abs(log10(EPSILON))) + 1)
    df.to_csv(fr'{WORKDIR}\Results\{method_name}.csv', index = False)

    print()
    print("Method:", method_name)
    print("n on last iteration:", n_last)

    
if __name__ == "__main__":
    EPSILON = 0.001
    a, b = 0.0, 0.5

    n_prev, x_prev, y_prev, n_last, x_last, y_last = calculate(euler_cauchy, f1)
    save("Euler_Cauchy_eq1", n_prev, x_prev, y_prev, n_last, x_last, y_last)

    n_prev, x_prev, y_prev, n_last, x_last, y_last = calculate(runge_kutta4, f1)
    save("Runge_Kutta4_eq1", n_prev, x_prev, y_prev, n_last, x_last, y_last)

    n_prev, x_prev, y_prev, n_last, x_last, y_last = calculate(adams3, f2_system)
    save("Adams3_eq2", n_prev, x_prev, y_prev, n_last, x_last, y_last)

    n_prev, x_prev, y_prev, n_last, x_last, y_last = calculate(adams4, f2_system)
    save("Adams4_eq2", n_prev, x_prev, y_prev, n_last, x_last, y_last)