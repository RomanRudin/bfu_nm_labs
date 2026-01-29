from lab import f, df, ddf, a, b
import numpy as np
print('\n'* 2)
print(f'Min(df(x)) on [a, b]: {(df(np.linspace(a, b, 10000))).min()}')
print(f'2/min(df(x)) on [a, b]: {2/(df(np.linspace(a, b, 10000))).min()}')