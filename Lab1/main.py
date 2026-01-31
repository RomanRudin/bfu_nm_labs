import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import *
from pathlib import Path
import os

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))

pd.options.display.float_format = '{:.8f}'.format
FUNCTION = 'exp(x - 1) + 2*x^2 - 7'                           # Change this to your function


def logger(file_name: str, data: list[int], interval: tuple[int, int]) -> pd.Series:
    result = pd.DataFrame({'x': data})
    with open(fr'{WORKDIR}\Results\{file_name}.csv', 'w') as file:
        result.to_csv(path_or_buf=file, sep='\t', lineterminator='\n')
    # print(result)
    return pd.DataFrame({'Метод решения': [file_name],
        'Выбранный интервал [a,b]': [str(interval)[1:-1]],
        'Полученное решение': [data[-1]],
        'Количество итераций': [len(data)]})

def f(x: float) -> float:
    return np.exp(x - 1) + 2*x**2 - 7   # Change this to your function

def df(x: float) -> float:
    return np.exp(x - 1) + 4*x          # Change this to your function's derivative

def ddf(x: float) -> float:
    return np.exp(x - 1) + 4            # Change this to your function's second derivative


a0, b0 = 0, 3                           # Change this to interval you want
n = (b0 - a0) * 100
step = 0.25                             # Change this to step you want
x_val = np.arange(a0, b0+step, step)
y_val = f(x_val)
print("Результаты расчетов")
print(f"a={a0}, b={b0}, n={(b0+step- a0) / step}")
print(pd.DataFrame({'x': x_val, 'f(x)': y_val}))
print("\n" * 2)
x_val = np.linspace(a0, b0, n)
y_val = f(x_val)


plt.figure()
plt.plot(x_val, y_val)
plt.axhline(0, color='black')
plt.grid()
plt.title(FUNCTION)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.savefig(r'Lab1/Results/1.png')


a, b = 1, 2                              # Change this to interval you want
eps = 1e-7                               # Change this to proximity you want
h = 1e-5                                 # Change this to h you want
full_result = pd.DataFrame(columns=['Метод решения', 'Выбранный интервал [a,b]', 'Полученное решение', 'Количество итераций'])


# Метод Ньютона
x = a if f(a) * ddf(a) > 0 else b
result = [f'{x:.8f}']
while abs(f(x)) > eps:
    x = x - f(x) / df(x)
    result.append(f'{x:.8f}')
full_result = pd.concat([full_result, logger('Newton method', result, (a, b))], ignore_index=True)


# Метод хорд
x = b
result = [f'{x:.8f}']
while abs(f(x)) > eps: 
    x = x - f(x) * (x - a) / (f(x) - f(a))
    result.append(f'{x:.8f}')
full_result = pd.concat([full_result, logger('Chords method', result, (a, b))], ignore_index=True)


# Метод секущих
x0, x1 = a, b
result = [f'{x0:.8f}']
while abs(f(x1)) > eps: 
    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0)) 
    x0, x1 = x1, x2
    result.append(f'{x1:.8f}') 
full_result = pd.concat([full_result, logger('Secants method', result, (a, b))], ignore_index=True)


# Конечноразностный метод Ньютона
x = a 
result = [f'{x:.8f}']
while abs(f(x)) > eps:
    x = x - h * f(x) / (f(x + h) - f(x))
    result.append(f'{x:.8f}')
full_result = pd.concat([full_result, logger('Finite differential Newton method', result, (a, b))], ignore_index=True)


# Метод Стеффенсена
x = 1.4                        # Change this to interval you want
result = [f'{x:.8f}']
while abs(f(x)) > eps: 
        x = x - f(x)**2 / (f(x + f(x)) - f(x))
        result.append(f'{x:.8f}')
full_result = pd.concat([full_result, logger('Steffensen method', result, (1.4, b))], ignore_index=True)


# Метод релаксации
theta = 0.1
x = b
max_iter = 10000 
result = [f'{x:.8f}']
for _ in range(max_iter): 
    x1 = x - theta * f(x)
    result.append(f'{x1:.8f}')
    if abs(x1 - x) < eps:
        break
    x = x1
full_result = pd.concat([full_result, logger('Relaxation method', result, (a, b))], ignore_index=True)


# Итоги
with open(fr'{WORKDIR}\Results\Results.csv', 'w', encoding='utf-8') as file:
    full_result.to_csv(path_or_buf=file, index=False, lineterminator='\n')
print(full_result)