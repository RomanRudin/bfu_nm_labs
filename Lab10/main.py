import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path


WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def f(x, y):
    return 6 * x - y


def initialize_grid(N):
    a, b = 0, 10
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = np.linspace(a, b, N + 1)
    U = np.zeros((N + 1, N + 1))
    for i in range(N + 1):
        U[i, 0] = x[i]
        U[i, N] = x[i] + 10
        U[0, i] = y[i]
        U[N, i] = 10 + y[i]
    return U, x, y, h



def jacobi(N, eps=0.01):
    U, x, y, h = initialize_grid(N)
    U_new = U.copy()
    diff = 1
    iterations = 0
    while diff > eps:
        diff = 0
        for i in range(1, N):
            for j in range(1, N):
                U_new[i, j] = 0.25 * (
                    U[i + 1, j] + U[i - 1, j] +
                    U[i, j + 1] + U[i, j - 1] -
                    h ** 2 * f(x[i], y[j])
                )
                diff = max(diff, abs(U_new[i, j] - U[i, j]))
        U[:] = U_new[:]
        iterations += 1
    return U, x, y, iterations


def zeidel(N, eps=0.01):
    U, x, y, h = initialize_grid(N)
    diff = 1
    iterations = 0
    while diff > eps:
        diff = 0
        for i in range(1, N):
            for j in range(1, N):
                old = U[i, j]
                U[i, j] = 0.25 * (
                    U[i + 1, j] + U[i - 1, j] +
                    U[i, j + 1] + U[i, j - 1] -
                    h ** 2 * f(x[i], y[j])
                )
                diff = max(diff, abs(U[i, j] - old))
        iterations += 1
    return U, x, y, iterations



def plot_3d(U, x, y, title):
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, U.T, cstride=1, edgecolor='none')
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('U')
    plt.savefig(rf"{WORKDIR}\Results\{title.replace(' ', '_')}.png", dpi=150, bbox_inches='tight')


def plot_slices(U, x, y, title):
    plt.figure()
    mid = len(y) // 2
    plt.plot(x, U[:, 0], label="y=0")
    plt.plot(x, U[:, mid], label=f"y={y[mid]:.1f}")
    plt.plot(x, U[:, -1], label="y=10")
    plt.legend()
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("U")
    plt.grid()
    plt.savefig(rf"{WORKDIR}\Results\{title.replace(' ', '_')}.png", dpi=150, bbox_inches='tight')



if __name__ == "__main__":
    grids = [5, 10]
    for N in grids:
        print()
        print(f"{N}x{N}")
        U_jacobi, x, y, it_j = jacobi(N)
        print(f"Jacobi method iterations: {it_j}")
        plot_3d(U_jacobi, x, y, f"3D Jacobi ({N}x{N})")
        plot_slices(U_jacobi, x, y, f"Jacobi slices ({N}x{N})")
        U_zeidel, x, y, it_s = zeidel(N)
        print(f"Zeidel method iterations: {it_s}")
        plot_3d(U_zeidel, x, y, f"3D Zeidel ({N}x{N})")
        plot_slices(U_zeidel, x, y, f"Zeidel slices ({N}x{N})")
