import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path


WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def phi(x):
    return x + 2 - x * (1 - x)


def g1(t) -> float:
    return 2.0


def g2(t) -> float:
    return 3.0



def solve_implicit_progonka(U_matrix, sigma):
    A_coeff = sigma * gamma
    B_coeff = sigma * gamma
    C_coeff = 1.0 + 2.0 * sigma * gamma

    for n in range(0, NT):
        alpha = np.zeros(NX + 1)
        beta = np.zeros(NX + 1)
        alpha[1] = 0.0
        beta[1] = g1(t_grid[n + 1])

        for i in range(1, NX):
            L_n = gamma * (U_matrix[n, i + 1] - 2 * U_matrix[n, i] + U_matrix[n, i - 1])
            F_i = U_matrix[n, i] + (1.0 - sigma) * L_n

            denom = C_coeff - A_coeff * alpha[i]
            alpha[i + 1] = B_coeff / denom
            beta[i + 1] = (F_i + A_coeff * beta[i]) / denom

        U_matrix[n + 1, NX] = g2(t_grid[n + 1])
        for i in range(NX - 1, 0, -1):
            U_matrix[n + 1, i] = alpha[i + 1] * U_matrix[n + 1, i + 1] + beta[i + 1]
        U_matrix[n + 1, 0] = g1(t_grid[n + 1])



def plot_2d(u, title):
    times_to_plot = [0, 2, 4, 6, 8, 10]
    indices = [int(np.argmin(np.abs(t_grid - t))) for t in times_to_plot]
    plt.figure(figsize=(10, 6))
    for t, idx in zip(times_to_plot, indices):
        plt.plot(x_grid, u[idx], label=f"t = {t}")
    plt.title(f"{title}")
    plt.xlabel("x")
    plt.ylabel("U(x,t)")
    plt.grid(True)
    plt.legend()
    plt.savefig(rf"{WORKDIR}\Results\{title.replace(' ', '_')}_2d.png", dpi=150, bbox_inches='tight')
    plt.close()


def plot_3d(u, title):
    X, T = np.meshgrid(x_grid, t_grid)
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    stride = int(NT / 100) if NT > 100 else 1
    surf = ax.plot_surface(X, T, u, rstride=stride, cstride=1, cmap='viridis', edgecolor='none')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('U(x,t)')
    plt.savefig(rf"{WORKDIR}\Results\{title.replace(' ', '_')}_3d.png", dpi=150, bbox_inches='tight')
    plt.close()



if __name__ == "__main__":
    A_bound, B_bound = 0.0, 1.0
    T0, TMAX = 0.0, 10.0
    D = 1.0
    NX = 20
    h = (B_bound - A_bound) / NX
    tau_max = (h ** 2) / (2.0 * D)
    NT = int(np.ceil((TMAX - T0) / tau_max)) + 50
    tau = (TMAX - T0) / NT
    x_grid = np.linspace(A_bound, B_bound, NX + 1)
    t_grid = np.linspace(T0, TMAX, NT + 1)
    gamma = D * tau / (h ** 2)

    U_expl = np.zeros((NT + 1, NX + 1))
    U_impl = np.zeros((NT + 1, NX + 1))

    for i in range(NX + 1):
        U_expl[0, i] = phi(x_grid[i])
        U_impl[0, i] = phi(x_grid[i])

    for n in range(0, NT):
        for i in range(1, NX):
            U_expl[n + 1, i] = U_expl[n, i] + gamma * (U_expl[n, i + 1] - 2 * U_expl[n, i] + U_expl[n, i - 1])
        U_expl[n + 1, 0] = g1(t_grid[n + 1])
        U_expl[n + 1, NX] = g2(t_grid[n + 1])

    solve_implicit_progonka(U_impl, 1.0)

    plot_2d(U_expl, "Явный метод")
    plot_2d(U_impl, "Неявный метод")
    plot_3d(U_expl, "Явный метод")
    plot_3d(U_impl, "Неявный метод")
    print("Графики сохранены в папку Results_9_lab")
