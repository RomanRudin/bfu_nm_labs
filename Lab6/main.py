import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
from pathlib import Path
from typing import Callable


WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def U0(x: float) -> float:
    return 2 * x**2 + 2

def f(x: float, t: float) -> float:
    return x

def U_left(t: float) -> float:
    return 2 * t**2 + 2

def U_right(t: float) -> float:
    return 2 * t**2 + 4



def runner(a: int, region: str, Nx: int, x: np.ndarray, idx_start: int, idx_end: int, bc_left: Callable, bc_right: Callable) -> list[bool | dict]:
    result = [False] * 4
    for index, scheme in enumerate(SCHEMES):
        U_sol = np.zeros((NT + 1, idx_end - idx_start + 1))
        U_old = U0(x)
        U_sol[0, :] = U_old[idx_start:idx_end + 1]
        t_curr = T0
        for j in range(NT):
            t_next = t_curr + dt
            U_new = scheme(a, region, U_old, Nx, x, bc_left, bc_right, t_curr, t_next)
            if U_new is False:
                result[index] = False
                break
            U_old = U_new
            t_curr = t_next
            U_sol[j+1, :] = U_new[idx_start:idx_end + 1]
        else:
            result[index] = {
                'idx_start' : idx_start, 
                'idx_end' : idx_end, 
                'U_sol' : U_sol
                }
    return result


def Scheme_1(a: int, region: str, U_old: np.ndarray, Nx: int, x: np.ndarray, bc_left: Callable, bc_right: Callable, t_curr: float, t_next: float) -> bool | np.ndarray[float]:
    if a < 0:
        return False
    U_new = np.zeros(Nx)
    if bc_left is not None:
        U_new[0] = bc_left(t_next)
    for i in range(1, Nx):
        U_new[i] = (U_old[i]
                    - (a * dt / dx) * (U_old[i] - U_old[i-1])
                    + dt * f(x[i], t_curr))
    return U_new

def Scheme_2(a: int, region: str, U_old: np.ndarray, Nx: int, x: np.ndarray, bc_left: Callable, bc_right: Callable, t_curr: float, t_next: float) -> bool | np.ndarray[float]:
    if a > 0:
        return False
    U_new = np.zeros(Nx)
    if bc_right is not None:
        U_new[Nx-1] = bc_right(t_next)
    for i in range(0, Nx-1):
        U_new[i] = (U_old[i]
                    - (a * dt / dx) * (U_old[i+1] - U_old[i])
                    + dt * f(x[i], t_curr))
    return U_new

def Scheme_3(a: int, region: str, U_old: np.ndarray, Nx: int, x: np.ndarray, bc_left: Callable, bc_right: Callable, t_curr: float, t_next: float) -> bool | np.ndarray[float]:
    if region == 'half' or a < 0:
        return False
    U_new = np.zeros(Nx)
    U_new[0] = bc_left(t_next)
    rhs = U_old + dt * f(x, t_next)
    for i in range(1, Nx):
        U_new[i] = (rhs[i] - ALPHA * U_new[i-1]) / BETA
    return U_new

def Scheme_4(a: int, region: str, U_old: np.ndarray, Nx: int, x: np.ndarray, bc_left: Callable, bc_right: Callable, t_curr: float, t_next: float) -> bool | np.ndarray[float]:
    if region == 'half':
        return False
    r = a * dt / dx
    U_new = np.zeros(Nx)
    U_new[0] = bc_left(t_next)
    t_mid = 0.5 * (t_curr + t_next)
    for i in range(1, Nx):
        x_mid = 0.5 * (x[i] + x[i - 1])
        F = dt * f(x_mid, t_mid)
        if abs(1 + r) < 1e-12:
            U_new[i] = U_old[i - 1] + F
        else:
            U_new[i] = ((r - 1) * U_new[i - 1]
                        + (1 - r) * U_old[i]
                        + (1 + r) * U_old[i - 1]
                        + 2 * F) / (1 + r)
    return U_new


def plot(a: int, region: str, results: list[bool | dict]) -> None:
    for scheme, result in enumerate(results):
        if result is False:
            continue
        idx_start, idx_end, U_sol = result['idx_start'], result['idx_end'], result['U_sol']
        X, T = np.meshgrid(x[idx_start:idx_end + 1], np.linspace(T0, T1, NT+1))
        fig = plt.figure(figsize = (10, 8))
        ax = fig.add_subplot(111, projection = '3d')
        ax.plot_surface(X, T, U_sol, cmap = 'viridis', edgecolor = 'none')
        ax.set_xlabel('x')
        ax.set_ylabel('t')
        ax.set_zlabel('U')
        ax.set_title(f'a = {a}, region = {region}, scheme = {scheme + 1}')
        filename = f"{WORKDIR}/Results/a{a}_{region}_{scheme + 1}.png"
        plt.savefig(filename, dpi = 150)
        plt.close(fig)



if __name__ == "__main__":
    L_extra = 30
    a_vals = [2, -2]
    dx = 0.1
    x0, x1 = 0.0, 1.0
    T0, T1 = 0.0, 10.0
    dt = 0.05
    NT = int((T1 - T0) / dt)

    SCHEMES = [Scheme_1, Scheme_2, Scheme_3, Scheme_4]


    for a, region in product(a_vals, ['half', 'rect']):
        ALPHA = -abs(a) * dt / dx
        BETA = 1 + abs(a) * dt / dx
        print(f"Computing: a={a}, region={region}")
        match region:
            case 'rect':
                x = np.linspace(x0, x1, int((x1 - x0) / dx) + 1)
                Nx = len(x)
                idx_start = 0
                idx_end = Nx - 1
                bc_left = U_left
                bc_right = U_right
            case 'half':
                if a > 0:
                    x_min = -L_extra
                    x_max = x1
                else:
                    x_min = x0
                    x_max = x1 + L_extra
                Nx = int(round((x_max - x_min) / dx)) + 1
                x = np.linspace(x_min, x_max, Nx)
                idx_start = int(round((0 - x_min) / dx))
                idx_end = int(round((1 - x_min) / dx))
                bc_left = lambda t, x0=x_min: U0(x0)
                bc_right = lambda t, x1=x_max: U0(x1)
        results = runner(a, region, Nx, x, idx_start, idx_end, bc_left, bc_right)
        plot(a, region, results)
    print(fr"Results are saved in {WORKDIR}\Results")