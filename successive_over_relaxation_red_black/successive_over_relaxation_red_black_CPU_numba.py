import numpy as np
from numba import njit, prange, get_num_threads
import sys
import time
# import matplotlib.pyplot as plt

# Uma iteração do método, compilada pelo Numba. Os vizinhos de um volume são
# sempre da outra cor, então os volumes de uma mesma cor não dependem uns dos
# outros e as linhas da malha podem ser divididas entre as threads
# (NUMBA_NUM_THREADS; sem a variável, todas)
# Os laços que andam de dois em dois volumes são escritos com while: com
# range(início, fim, 2) o código compilado pelo Numba fica mais lento
@njit(parallel=True)
def iteracao(Ap, Aw, Ae, As, An, Bp, phi_new, norma_b, omega):
    nVolY, nVolX = phi_new.shape

    # Volumes vermelhos (i + j par) e depois pretos (i + j ímpar)
    for cor in range(2):
        # Face oeste
        j = 0
        i = 2 - cor
        while i < nVolY - 1:
            phi_new[i, j] = (1 - omega) * phi_new[i, j] + omega * (
                - Ae[i, j] * phi_new[i, j + 1]
                - As[i, j] * phi_new[i - 1, j]
                - An[i, j] * phi_new[i + 1, j]
                + Bp[i, j]
            ) / Ap[i, j]
            i += 2

        # Face leste
        j = nVolX - 1
        i = 2 - (cor + nVolX - 1) % 2
        while i < nVolY - 1:
            phi_new[i, j] = (1 - omega) * phi_new[i, j] + omega * (
                - Aw[i, j] * phi_new[i, j - 1]
                - As[i, j] * phi_new[i - 1, j]
                - An[i, j] * phi_new[i + 1, j]
                + Bp[i, j]
            ) / Ap[i, j]
            i += 2

        # Face sul
        i = 0
        j = 2 - cor
        while j < nVolX - 1:
            phi_new[i, j] = (1 - omega) * phi_new[i, j] + omega * (
                - Aw[i, j] * phi_new[i, j - 1]
                - Ae[i, j] * phi_new[i, j + 1]
                - An[i, j] * phi_new[i + 1, j]
                + Bp[i, j]
            ) / Ap[i, j]
            j += 2

        # Face norte
        i = nVolY - 1
        j = 2 - (cor + nVolY - 1) % 2
        while j < nVolX - 1:
            phi_new[i, j] = (1 - omega) * phi_new[i, j] + omega * (
                - Aw[i, j] * phi_new[i, j - 1]
                - Ae[i, j] * phi_new[i, j + 1]
                - As[i, j] * phi_new[i - 1, j]
                + Bp[i, j]
            ) / Ap[i, j]
            j += 2

        # Volumes internos: em cada linha, só os volumes da cor atual
        for i in prange(1, nVolY - 1):
            j = 2 - (cor + i) % 2
            while j < nVolX - 1:
                phi_new[i, j] = (1 - omega) * phi_new[i, j] + omega * (
                    - Aw[i, j] * phi_new[i, j - 1]
                    - Ae[i, j] * phi_new[i, j + 1]
                    - As[i, j] * phi_new[i - 1, j]
                    - An[i, j] * phi_new[i + 1, j]
                    + Bp[i, j]
                ) / Ap[i, j]
                j += 2

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

# Centros dos volumes na direção y. A solução não depende da posição em x
y = np.linspace(-R - dx / 2, R + dx / 2, nVolY + 2)

nVolX = int(np.round(L / dx))
L = nVolX * dx

nVolX += 2
nVolY += 2

# Campo de velocidade: só depende de y, e é o mesmo nas faces Oeste e Leste e
# no centro do volume. Não há velocidade na direção y
dPdx = (Pout - Pin) / L

u = 1 / (4 * mu) * (-dPdx) * (R**2 - y**2)

vs = 0.0
vn = 0.0

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

# Volumes internos: a velocidade de cada linha vale para todas as colunas
ui = u[1:-1, None]

Ap[1:-1, 1:-1] = dx * rho * (- np.minimum(0, ui) + np.maximum(0, ui) - min(0, vs) + max(0, vn)) + 4 * Gamma

Aw[1:-1, 1:-1] = -dx * rho * np.maximum(0, ui) - Gamma
Ae[1:-1, 1:-1] = dx * rho * np.minimum(0, ui) - Gamma
As[1:-1, 1:-1] = -dx * rho * max(0, vs) - Gamma
An[1:-1, 1:-1] = dx * rho * min(0, vn) - Gamma
Bp[1:-1, 1:-1] = 0

# Solução inicial de phi
phi_new = np.zeros((nVolY, nVolX))

# Resolução do sistema linear
residuo_iteracao = 1
numero_iteracao = 0
numero_maximo_iteracao = 1000000
residuo_final = 1e-10
omega = 1.0  # Fator de relaxação (1 = Gauss-Seidel red-black)

# Norma do termo fonte (os cantos fantasmas não fazem parte do sistema)
b = Bp.copy()
b[::nVolY - 1, ::nVolX - 1] = 0
norma_b = np.linalg.norm(b)

print("=> Início das iterações")

# A primeira chamada da função compila o código; ela é descartada
aquecimento = True

end_allocation_time = time.perf_counter()
start_iteration_time = end_allocation_time

while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
    residuo_iteracao = iteracao(Ap, Aw, Ae, As, An, Bp, phi_new, norma_b, omega)
    numero_iteracao += 1

    if aquecimento:
        aquecimento = False
        phi_new[:] = 0
        residuo_iteracao = 1
        numero_iteracao = 0
        end_warm_up_time = time.perf_counter()
        start_iteration_time = end_warm_up_time

end_time = time.perf_counter()

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
# x = np.linspace(0 - dx / 2, L + dx / 2, nVolX)
# X, Y = np.meshgrid(x, y)
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
