import numpy as np
from numba import njit, prange, get_num_threads
import sys
import time
# import matplotlib.pyplot as plt

# Uma iteração do método, compilada pelo Numba. As linhas da malha são
# divididas entre as threads (NUMBA_NUM_THREADS; sem a variável, todas)
@njit(parallel=True)
def iteracao(Ap, Aw, Ae, As, An, Bp, phi_old, phi_new, norma_b):
    nVolY, nVolX = phi_old.shape

    # Atualiza as bordas Oeste
    j = 0
    for i in range(1, nVolY - 1):
        phi_old[i, j] = (
            - Ae[i, j] * phi_old[i, j + 1]
            - As[i, j] * phi_old[i - 1, j]
            - An[i, j] * phi_old[i + 1, j]
            + Bp[i, j]
        ) / Ap[i, j]

    # Atualiza as bordas Leste
    j = nVolX - 1
    for i in range(1, nVolY - 1):
        phi_old[i, j] = (
            - Aw[i, j] * phi_old[i, j - 1]
            - As[i, j] * phi_old[i - 1, j]
            - An[i, j] * phi_old[i + 1, j]
            + Bp[i, j]
        ) / Ap[i, j]

    # Atualiza as bordas Sul
    i = 0
    for j in range(1, nVolX - 1):
        phi_old[i, j] = (
            - Aw[i, j] * phi_old[i, j - 1]
            - Ae[i, j] * phi_old[i, j + 1]
            - An[i, j] * phi_old[i + 1, j]
            + Bp[i, j]
        ) / Ap[i, j]

    # Atualiza as bordas Norte
    i = nVolY - 1
    for j in range(1, nVolX - 1):
        phi_old[i, j] = (
            - Aw[i, j] * phi_old[i, j - 1]
            - Ae[i, j] * phi_old[i, j + 1]
            - As[i, j] * phi_old[i - 1, j]
            + Bp[i, j]
        ) / Ap[i, j]

    # As bordas (e os cantos) valem para os dois vetores
    for i in range(nVolY):
        phi_new[i, 0] = phi_old[i, 0]
        phi_new[i, nVolX - 1] = phi_old[i, nVolX - 1]

    for j in range(nVolX):
        phi_new[0, j] = phi_old[0, j]
        phi_new[nVolY - 1, j] = phi_old[nVolY - 1, j]

    # Volumes internos: lê só a iteração anterior (phi_old), o que mantém o
    # método como Jacobi. Com um vetor só, o laço seria Gauss-Seidel
    for i in prange(1, nVolY - 1):
        for j in range(1, nVolX - 1):
            phi_new[i, j] = (- Aw[i, j] * phi_old[i, j - 1]
                - Ae[i, j] * phi_old[i, j + 1]
                - As[i, j] * phi_old[i - 1, j]
                - An[i, j] * phi_old[i + 1, j]
                + Bp[i, j]) / Ap[i, j]

    # Resíduo real do sistema linear: ||b - A*phi|| / ||b||
    # (os cantos fantasmas não fazem parte do sistema)
    soma = 0.0

    for i in prange(nVolY):
        for j in range(nVolX):
            if (i == 0 or i == nVolY - 1) and (j == 0 or j == nVolX - 1):
                continue

            res = Bp[i, j] - Ap[i, j] * phi_new[i, j]

            if j > 0:
                res -= Aw[i, j] * phi_new[i, j - 1]

            if j < nVolX - 1:
                res -= Ae[i, j] * phi_new[i, j + 1]

            if i > 0:
                res -= As[i, j] * phi_new[i - 1, j]

            if i < nVolY - 1:
                res -= An[i, j] * phi_new[i + 1, j]

            soma += res * res

    return np.sqrt(soma) / norma_b

# Início da medição do tempo de alocação
start_allocation_time = time.perf_counter()

# Configuração inicial
L = 2.0  # Comprimento do domínio [m]
D = 0.01  # Diâmetro do canal [m]
rho = 1e3  # Densidade da água [kg/m³]
mu = 8.9e-4  # Viscosidade da água [kg/m·s]
Gamma = 0.61 / 4200  # Coeficiente de condução [W/(m·K)]
Pin = 1.0  # Pressão na entrada [Pa]
Pout = 0.0  # Pressão na saída [Pa]
Tin = 25.0  # Temperatura de entrada [°C]
Twall = 100.0  # Temperatura nas paredes [°C]

# Número de volumes na direção y (primeiro argumento da linha de comando)
nVolY = int(sys.argv[1]) if len(sys.argv) > 1 else 400

# Configuração da malha
R = D / 2
dx = D / nVolY
y = np.linspace(-R - dx / 2, R + dx / 2, nVolY + 2)

nVolX = int(np.round(L / dx))
L = nVolX * dx
x = np.linspace(0 - dx / 2, L + dx / 2, nVolX + 2)

nVolX += 2
nVolY += 2

X, Y = np.meshgrid(x, y)

# Campo de velocidade
dPdx = (Pout - Pin) / L

xw = X - dx / 2
yw = Y
uw = 1 / (4 * mu) * (-dPdx) * (R**2 - yw**2)

xe = X + dx / 2
ye = Y
ue = 1 / (4 * mu) * (-dPdx) * (R**2 - ye**2)

vs = np.zeros_like(Y)
vn = np.zeros_like(Y)

up = 1 / (4 * mu) * (-dPdx) * (R**2 - Y**2)
vp = np.zeros_like(Y)
Vp = np.sqrt(up**2 + vp**2)

# Coeficientes da equação de temperatura
Ap = np.zeros((nVolY, nVolX))
Aw = np.zeros((nVolY, nVolX))
Ae = np.zeros((nVolY, nVolX))
As = np.zeros((nVolY, nVolX))
An = np.zeros((nVolY, nVolX))
Bp = np.zeros((nVolY, nVolX))

# Faces Oeste
Ap[:, 0] = 1
Aw[:, 0] = 0
Ae[:, 0] = 1
As[:, 0] = 0
An[:, 0] = 0
Bp[:, 0] = 2 * Tin

# Faces Leste
Ap[:, -1] = 1
Aw[:, -1] = -1
Ae[:, -1] = 0
As[:, -1] = 0
An[:, -1] = 0
Bp[:, -1] = 0

# Faces Sul
Ap[0, :] = 1
Aw[0, :] = 0
Ae[0, :] = 0
As[0, :] = 0
An[0, :] = 1
Bp[0, :] = 2 * Twall

# Faces Norte
Ap[-1, :] = 1
Aw[-1, :] = 0
Ae[-1, :] = 0
As[-1, :] = 1
An[-1, :] = 0
Bp[-1, :] = 2 * Twall

# Volumes internos
Ap[1:-1, 1:-1] = dx * rho * (- np.minimum(0, uw[1:-1, 1:-1]) + np.maximum(0, ue[1:-1, 1:-1]) - np.minimum(0, vs[1:-1, 1:-1]) + np.maximum(0, vn[1:-1, 1:-1])) + 4 * Gamma

Aw[1:-1, 1:-1] = -dx * rho * np.maximum(0, uw[1:-1, 1:-1]) - Gamma
Ae[1:-1, 1:-1] = dx * rho * np.minimum(0, ue[1:-1, 1:-1]) - Gamma
As[1:-1, 1:-1] = -dx * rho * np.maximum(0, vs[1:-1, 1:-1]) - Gamma
An[1:-1, 1:-1] = dx * rho * np.minimum(0, vn[1:-1, 1:-1]) - Gamma
Bp[1:-1, 1:-1] = 0

# Solução inicial de phi
phi_new = np.zeros((nVolY, nVolX))

# Resolução do sistema linear
residuo_iteracao = 1
numero_iteracao = 0
numero_maximo_iteracao = 1000000
residuo_final = 1e-10

# Norma do termo fonte (os cantos fantasmas não fazem parte do sistema)
b = Bp.copy()
b[::nVolY - 1, ::nVolX - 1] = 0
norma_b = np.linalg.norm(b)

# O laço atualiza um vetor a partir do outro, então a solução inicial é
# necessária nos dois
phi_old = phi_new.copy()

print("=> Início das iterações")

# A primeira chamada da função compila o código; ela é descartada
aquecimento = True

end_allocation_time = time.perf_counter()
start_iteration_time = end_allocation_time

while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
    residuo_iteracao = iteracao(Ap, Aw, Ae, As, An, Bp, phi_old, phi_new, norma_b)
    phi_old, phi_new = phi_new, phi_old
    numero_iteracao += 1

    if aquecimento:
        aquecimento = False
        phi_old[:] = 0
        phi_new[:] = 0
        residuo_iteracao = 1
        numero_iteracao = 0
        end_warm_up_time = time.perf_counter()
        start_iteration_time = end_warm_up_time

end_time = time.perf_counter()

# Depois da troca, a solução da última iteração está em phi_old
phi_new = phi_old

estado = "convergiu" if residuo_iteracao <= residuo_final else "nao_convergiu"

print(f"=> Resultado: nVolY={nVolY - 2}"
      f" alocacao={end_allocation_time - start_allocation_time:.6f}"
      f" iteracao={end_time - start_iteration_time:.6f}"
      f" numero_iteracao={numero_iteracao}"
      f" residuo={float(residuo_iteracao):.6e}"
      f" estado={estado}"
      f" aquecimento={end_warm_up_time - end_allocation_time:.6f}"
      f" threads={get_num_threads()}")

# Exibição dos resultados
# plt.figure()
# plt.contourf(
#     X[1:-1, 1:-1],
#     Y[1:-1, 1:-1],
#     phi_new[1:-1, 1:-1],
#     cmap="jet",
# )
# plt.colorbar(label="Temperatura (°C)")
# plt.title("Campo de Temperatura")
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# plt.axis("equal")
# plt.savefig("campo_temperatura.png")
# print("=> Gráfico salvo como 'campo_temperatura.png'")
