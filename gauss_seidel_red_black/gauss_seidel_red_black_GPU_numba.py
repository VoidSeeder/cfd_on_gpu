import numpy as np
from numba_cuda_mlir import cuda
import sys
import time
# import matplotlib.pyplot as plt

# Kernels de uma iteração do método, compilados para a GPU. Os vizinhos de um
# volume são sempre da outra cor, então os volumes de uma mesma cor não
# dependem uns dos outros: cada thread da GPU calcula um volume da cor atual
# (ou, no resíduo, uma coluna) e grava no próprio phi_new

# Faces da cor atual (0 = vermelho, i + j par; 1 = preto, i + j ímpar). A
# thread k calcula as faces Oeste e Leste da linha k e as faces Sul e Norte da
# coluna k
@cuda.jit
def bordas(Ap, Aw, Ae, As, An, Bp, phi_new, cor):
    nVolY, nVolX = phi_new.shape
    k = cuda.grid(1)

    if k >= 1 and k < nVolY - 1:
        # Atualiza as bordas Oeste
        i = k
        j = 0
        if (i + j) % 2 == cor:
            phi_new[i, j] = (
                - Ae[i, j] * phi_new[i, j + 1]
                - As[i, j] * phi_new[i - 1, j]
                - An[i, j] * phi_new[i + 1, j]
                + Bp[i, j]
            ) / Ap[i, j]

        # Atualiza as bordas Leste
        j = nVolX - 1
        if (i + j) % 2 == cor:
            phi_new[i, j] = (
                - Aw[i, j] * phi_new[i, j - 1]
                - As[i, j] * phi_new[i - 1, j]
                - An[i, j] * phi_new[i + 1, j]
                + Bp[i, j]
            ) / Ap[i, j]

    if k >= 1 and k < nVolX - 1:
        # Atualiza as bordas Sul
        i = 0
        j = k
        if (i + j) % 2 == cor:
            phi_new[i, j] = (
                - Aw[i, j] * phi_new[i, j - 1]
                - Ae[i, j] * phi_new[i, j + 1]
                - An[i, j] * phi_new[i + 1, j]
                + Bp[i, j]
            ) / Ap[i, j]

        # Atualiza as bordas Norte
        i = nVolY - 1
        if (i + j) % 2 == cor:
            phi_new[i, j] = (
                - Aw[i, j] * phi_new[i, j - 1]
                - Ae[i, j] * phi_new[i, j + 1]
                - As[i, j] * phi_new[i - 1, j]
                + Bp[i, j]
            ) / Ap[i, j]

# Volumes internos da cor atual: a thread (k, i) calcula o k-ésimo volume
# dessa cor na linha i, então nenhuma thread fica sem volume para calcular
@cuda.jit
def internos(Ap, Aw, Ae, As, An, Bp, phi_new, cor):
    nVolY, nVolX = phi_new.shape
    k, i = cuda.grid(2)
    j = 2 - (cor + i) % 2 + 2 * k

    if i >= 1 and i < nVolY - 1 and j < nVolX - 1:
        phi_new[i, j] = (
            - Aw[i, j] * phi_new[i, j - 1]
            - Ae[i, j] * phi_new[i, j + 1]
            - As[i, j] * phi_new[i - 1, j]
            - An[i, j] * phi_new[i + 1, j]
            + Bp[i, j]
        ) / Ap[i, j]

# Resíduo real do sistema linear, b - A*phi: a thread j soma os quadrados do
# resíduo dos volumes da coluna j (os cantos fantasmas não fazem parte do
# sistema)
@cuda.jit
def residuo(Ap, Aw, Ae, As, An, Bp, phi_new, soma_coluna):
    nVolY, nVolX = phi_new.shape
    j = cuda.grid(1)

    if j < nVolX:
        soma = 0.0

        for i in range(nVolY):
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

        soma_coluna[j] = soma

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

# Norma do termo fonte (os cantos fantasmas não fazem parte do sistema)
b = Bp.copy()
b[::nVolY - 1, ::nVolX - 1] = 0
norma_b = np.linalg.norm(b)

# Os coeficientes são montados na memória principal e copiados para a GPU
Ap_gpu = cuda.to_device(Ap)
Aw_gpu = cuda.to_device(Aw)
Ae_gpu = cuda.to_device(Ae)
As_gpu = cuda.to_device(As)
An_gpu = cuda.to_device(An)
Bp_gpu = cuda.to_device(Bp)

# As duas cores são atualizadas no mesmo vetor
phi_new_gpu = cuda.to_device(phi_new)

# Soma dos quadrados do resíduo de cada coluna, na GPU e na memória principal
soma_coluna_gpu = cuda.device_array(nVolX)
soma_coluna = np.zeros(nVolX)

# Divisão das threads da GPU em blocos: uma dimensão para as bordas e o
# resíduo, duas para os volumes internos (em cada linha, só metade dos volumes
# é de cada cor)
threads_1d = 256
blocos_1d = (max(nVolX, nVolY) + threads_1d - 1) // threads_1d

threads_2d = (16, 16)
blocos_2d = ((nVolX // 2 + threads_2d[0] - 1) // threads_2d[0], (nVolY + threads_2d[1] - 1) // threads_2d[1])

print("=> Início das iterações")

# A primeira iteração compila os kernels; ela é descartada
aquecimento = True

cuda.synchronize()
end_allocation_time = time.perf_counter()
start_iteration_time = end_allocation_time

while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
    # Volumes vermelhos (i + j par) e depois pretos (i + j ímpar)
    for cor in range(2):
        bordas[blocos_1d, threads_1d](Ap_gpu, Aw_gpu, Ae_gpu, As_gpu, An_gpu, Bp_gpu, phi_new_gpu, cor)
        internos[blocos_2d, threads_2d](Ap_gpu, Aw_gpu, Ae_gpu, As_gpu, An_gpu, Bp_gpu, phi_new_gpu, cor)

    residuo[blocos_1d, threads_1d](Ap_gpu, Aw_gpu, Ae_gpu, As_gpu, An_gpu, Bp_gpu, phi_new_gpu, soma_coluna_gpu)

    # ||b - A*phi|| / ||b||: as somas das colunas são copiadas da GPU e
    # somadas na CPU
    soma_coluna_gpu.copy_to_host(soma_coluna)
    residuo_iteracao = np.sqrt(np.sum(soma_coluna)) / norma_b

    numero_iteracao += 1

    if aquecimento:
        aquecimento = False
        phi_new_gpu.copy_to_device(phi_new)
        residuo_iteracao = 1
        numero_iteracao = 0
        cuda.synchronize()
        end_warm_up_time = time.perf_counter()
        start_iteration_time = end_warm_up_time

cuda.synchronize()
end_time = time.perf_counter()

phi_new = phi_new_gpu.copy_to_host()

estado = "convergiu" if residuo_iteracao <= residuo_final else "nao_convergiu"

print(f"=> Resultado: nVolY={nVolY - 2}"
      f" alocacao={end_allocation_time - start_allocation_time:.6f}"
      f" iteracao={end_time - start_iteration_time:.6f}"
      f" numero_iteracao={numero_iteracao}"
      f" residuo={float(residuo_iteracao):.6e}"
      f" estado={estado}"
      f" aquecimento={end_warm_up_time - end_allocation_time:.6f}")

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
# plt.savefig("campo_temperatura_GPU.png")
# print("=> Gráfico salvo como 'campo_temperatura.png'")
