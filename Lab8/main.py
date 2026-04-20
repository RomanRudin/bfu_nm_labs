import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

WORKDIR = Path(__file__).parent
if not Path.exists(Path(fr'{WORKDIR}\Results')):
    Path.mkdir(Path(fr'{WORKDIR}\Results'))


def phi(x):   return 1.0 - x
def psi(x):   return -1.0 
def g1(t):    return 1.0  
def g2(t):    return 0.0  



def _init() -> np.ndarray:
    U = np.zeros((NT + 1, NX + 1))
    
    U[0, :] = phi(x)
    U[0, 0], U[0, NX] = g1(t[0]), g2(t[0])

    U[1, 1:-1] = U[0, 1:-1] + tau * psi(x[1:-1]) 
    U[1, 0],  U[1, NX] = g1(t[1]), g2(t[1])
    return U



def solve_explicit() -> np.ndarray:
    U = _init()
    for j in range(1, NT):
        i = np.arange(1, NX)                         
        U[j+1, i] = (2*U[j, i] - U[j-1, i]
                     + r2 * (U[j, i+1] - 2*U[j, i] + U[j, i-1]))
        U[j+1, 0]  = g1(t[j+1])
        U[j+1, NX] = g2(t[j+1])
    return U 



def progonka(F: np.ndarray, A: float, B: float, C: float, bc_left: float, bc_right: float) -> np.ndarray:
    alpha = np.zeros(NX + 1)
    beta = np.zeros(NX + 1)

    alpha[1] = 0.0
    beta[1] = bc_left
    for i in range(1, NX):          
        denom = C - A * alpha[i]
        alpha[i+1] = B / denom
        beta[i+1] = (F[i-1] + A * beta[i]) / denom

    u = np.empty(NX + 1)
    u[NX] = bc_right
    for i in range(NX - 1, 0, -1):
        u[i] = alpha[i+1] * u[i+1] + beta[i+1]
    u[0] = bc_left
    return u



def solve_implicit_T() -> np.ndarray:
    A = r2
    C = 1.0 + 2.0 * r2
    U = _init()
    for j in range(1, NT):
        F = 2*U[j, 1:-1] - U[j-1, 1:-1]

        bl = g1(t[j+1])
        br = g2(t[j+1])
        F[0]  += A * bl
        F[-1] += A * br

        U[j+1, :] = progonka(F, A, A, C, bl, br)
    return U



def solve_implicit_2() -> np.ndarray:
    sigma = 0.25
    A = sigma * r2
    C = 1.0 + 2.0 * sigma * r2
    U = _init()

    for j in range(1, NT):
        Lap_j  = U[j, 2:] - 2*U[j, 1:-1] + U[j, :-2]
        Lap_jm = U[j-1, 2:] - 2*U[j-1, 1:-1] + U[j-1, :-2]

        F = (2*U[j, 1:-1] - U[j-1, 1:-1] + (1 - 2*sigma) * r2 * Lap_j + sigma * r2 * Lap_jm)
        bl = g1(t[j+1])
        br = g2(t[j+1])
        F[0] += A * bl
        F[-1] += A * br

        U[j+1, :] = progonka(F, A, A, C, bl, br)
    return U



def plot_results(u, title) -> None:
    t = np.linspace(0.0, 10.0, NT + 1)
    times_to_plot = [0, 2, 4, 6, 8, 10]
    indices = [int(np.argmin(np.abs(t - tc))) for tc in times_to_plot]
    plt.figure(figsize = (10, 6))
    for t, idx in zip(times_to_plot, indices):
        plt.plot(x, u[idx], label=f"t = {t}")
    plt.title(f"{title} (2D)")
    plt.xlabel("x")
    plt.ylabel("U(x,t)")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"{WORKDIR}/Results/{title.replace(' ', '_')}_2d.png")
    plt.close()
    t = np.linspace(0.0, 10.0, NT + 1)
    X, T = np.meshgrid(x, t)
    fig = plt.figure(figsize = (12, 8))
    ax = fig.add_subplot(111, projection = '3d')
    surf = ax.plot_surface(X, T, u, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_xlabel("x"); ax.set_ylabel("t"); ax.set_zlabel("U")
    fig.colorbar(surf, ax = ax, shrink = 0.5, aspect = 5)
    ax.set_title(f"{title} (3D)")
    plt.savefig(f"{WORKDIR}/Results/{title.replace(' ', '_')}_3d.png")
    plt.close()


if __name__ == "__main__":
    D = 1.0
    NX = 20                       
    h = 1.0 / NX                 
    tau = h / D                    
    NT = int(np.ceil(10.0 / tau)) 
    tau = 10.0 / NT                
    r2 = (D * tau / h) ** 2       

    x = np.linspace(0.0, 1.0, NX + 1)
    t = np.linspace(0.0, 10.0, NT + 1)


    U_expl = solve_explicit()
    U_impl1 = solve_implicit_T()
    U_impl2 = solve_implicit_2()



    plot_results(U_expl, "Явный метод")
    plot_results(U_impl1, "Неявный метод 1")
    plot_results(U_impl2, "Неявный метод 2")

    print(f"Графики сохранены в папку {WORKDIR}/Results/")
