import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))



def initial_condition(x) -> np.ndarray:
    return np.where(x < 0.5, 4.0, 2.0)

def exact_solution(x: np.ndarray, t: float) -> np.ndarray:
    shock_pos = 0.5 + 3.0 * t
    if shock_pos <= 1.0:
        return np.where(x < shock_pos, 4.0, 2.0)
    return np.full_like(x, 4.0)

def flux(u) -> float:
    return 0.5 * u ** 2


def solve_viscosity() -> np.ndarray:
    U = initial_condition(x)
    nu = 0.01
    for _ in range(nt):
        U_new = U.copy()
        for i in range(1, nx - 1):
            convection = (flux(U[i + 1]) - flux(U[i - 1])) / (2 * dx)
            diffusion = nu * (U[i + 1] - 2 * U[i] + U[i - 1]) / dx ** 2
            U_new[i] = U[i] - dt * convection + dt * diffusion
        U_new[0] = U_new[1]
        U_new[-1] = U_new[-2]
        U = U_new.copy()
    return U


def solve_conservative() -> np.ndarray:
    U = initial_condition(x)
    for _ in range(nt):
        U_new = U.copy()
        for i in range(1, nx - 1):
            U_new[i] = 0.5 * (U[i + 1] + U[i - 1]) \
                       - dt / (2 * dx) * (flux(U[i + 1]) - flux(U[i - 1]))
        U_new[0] = U_new[1]
        U_new[-1] = U_new[-2]
        U = U_new.copy()
    return U


def solve_viscosity_3d() -> np.ndarray:
    U = initial_condition(x)
    nu = 0.01
    history = [U.copy()]
    for _ in range(nt):
        U_new = U.copy()
        for i in range(1, nx - 1):
            convection = (flux(U[i + 1]) - flux(U[i - 1])) / (2 * dx)
            diffusion = nu * (U[i + 1] - 2 * U[i] + U[i - 1]) / dx ** 2
            U_new[i] = U[i] - dt * convection + dt * diffusion
        U_new[0] = U_new[1]
        U_new[-1] = U_new[-2]
        U = U_new.copy()
        history.append(U.copy())
    return np.array(history)

def solve_conservative_3d() -> np.ndarray:
    U = initial_condition(x)
    history = [U.copy()]
    for _ in range(nt):
        U_new = U.copy()
        for i in range(1, nx - 1):
            U_new[i] = 0.5 * (U[i + 1] + U[i - 1]) \
                       - dt / (2 * dx) * (flux(U[i + 1]) - flux(U[i - 1]))
        U_new[0] = U_new[1]
        U_new[-1] = U_new[-2]
        U = U_new.copy()
        history.append(U.copy())
    return np.array(history)


def plot(U_init: np.ndarray, U_con: np.ndarray, u_exact: np.ndarray, title: str, filename: str) -> None:
    plt.figure(figsize=(8, 5))
    plt.plot(x, U_init, '--', label="Начальное условие")
    plt.plot(x, U_con, label=title)
    plt.plot(x, u_exact, 'b-', lw=2, label='Точное решение')
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("U")
    plt.grid()
    plt.legend()
    plt.savefig(f"{WORKDIR}/Results/{filename}.png", dpi=150)


def plot3d(U_3d: np.ndarray, title: str, filename: str) -> None:
    t = np.linspace(0, nt * dt, nt + 1)
    X, T = np.meshgrid(x, t)
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, U_3d, cmap='viridis', edgecolor='none', alpha=0.9)
    ax.set_title(f"{title} (3D)")
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_zlabel("U")
    plt.savefig(f"{WORKDIR}/Results/{filename}3d.png", dpi=150)




if __name__ == "__main__":
    nx = 200
    x_start, x_end = 0, 1
    x = np.linspace(x_start, x_end, nx)
    dx = (x_end - x_start) / nx
    Umax = 4
    dt = 0.4 * dx / Umax  # CFL
    nt = 300

    U_visc = solve_viscosity()
    U_con = solve_conservative()
    U_init = initial_condition(x)
    U_exact = exact_solution(x, 0.15)
    U_visc_3d = solve_viscosity_3d()
    U_con_3d = solve_conservative_3d()

    plot(U_init, U_con, U_exact, "Консервативная схема", "conservative")
    plot(U_init, U_visc, U_exact, "Искусственная вязкость", "viscosity")
    plot3d(U_con_3d, "Консервативная схема", "conservative_3d")
    plot3d(U_visc_3d, "Искусственная вязкость", "viscosity_3d")