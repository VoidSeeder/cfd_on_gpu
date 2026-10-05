import numpy as np
import sys
import time
# import matplotlib.pyplot as plt

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

print("=> Início das iterações")

end_allocation_time = time.perf_counter()
start_iteration_time = end_allocation_time

while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
    # Volumes vermelhos (i + j par) e depois pretos (i + j ímpar).
    # Os vizinhos de um volume são sempre da outra cor.
    for cor in (0, 1):
        # Face oeste
        i = 2 - cor
        phi_new[i:-1:2, 0] = (
            - Ae[i:-1:2, 0] * phi_new[i:-1:2, 1]
            - As[i:-1:2, 0] * phi_new[i - 1:-2:2, 0]
            - An[i:-1:2, 0] * phi_new[i + 1::2, 0]
            + Bp[i:-1:2, 0]
        ) / Ap[i:-1:2, 0]

        # Face leste
        i = 2 - (cor + nVolX - 1) % 2
        phi_new[i:-1:2, -1] = (
            - Aw[i:-1:2, -1] * phi_new[i:-1:2, -2]
            - As[i:-1:2, -1] * phi_new[i - 1:-2:2, -1]
            - An[i:-1:2, -1] * phi_new[i + 1::2, -1]
            + Bp[i:-1:2, -1]
        ) / Ap[i:-1:2, -1]

        # Face sul
        j = 2 - cor
        phi_new[0, j:-1:2] = (
            - Aw[0, j:-1:2] * phi_new[0, j - 1:-2:2]
            - Ae[0, j:-1:2] * phi_new[0, j + 1::2]
            - An[0, j:-1:2] * phi_new[1, j:-1:2]
            + Bp[0, j:-1:2]
        ) / Ap[0, j:-1:2]

        # Face norte
        j = 2 - (cor + nVolY - 1) % 2
        phi_new[-1, j:-1:2] = (
            - Aw[-1, j:-1:2] * phi_new[-1, j - 1:-2:2]
            - Ae[-1, j:-1:2] * phi_new[-1, j + 1::2]
            - As[-1, j:-1:2] * phi_new[-2, j:-1:2]
            + Bp[-1, j:-1:2]
        ) / Ap[-1, j:-1:2]

        # Volumes internos: linhas ímpares e depois linhas pares
        for i in (1, 2):
            j = 2 - (cor + i) % 2
            phi_new[i:-1:2, j:-1:2] = (
                - Aw[i:-1:2, j:-1:2] * phi_new[i:-1:2, j - 1:-2:2]
                - Ae[i:-1:2, j:-1:2] * phi_new[i:-1:2, j + 1::2]
                - As[i:-1:2, j:-1:2] * phi_new[i - 1:-2:2, j:-1:2]
                - An[i:-1:2, j:-1:2] * phi_new[i + 1::2, j:-1:2]
                + Bp[i:-1:2, j:-1:2]
            ) / Ap[i:-1:2, j:-1:2]

    # Resíduo real do sistema linear: ||b - A*phi|| / ||b||
    res = Bp - Ap * phi_new
    res[:, 1:] -= Aw[:, 1:] * phi_new[:, :-1]
    res[:, :-1] -= Ae[:, :-1] * phi_new[:, 1:]
    res[1:, :] -= As[1:, :] * phi_new[:-1, :]
    res[:-1, :] -= An[:-1, :] * phi_new[1:, :]
    res[::nVolY - 1, ::nVolX - 1] = 0

    residuo_iteracao = np.linalg.norm(res) / norma_b
    numero_iteracao += 1

end_time = time.perf_counter()

estado = "convergiu" if residuo_iteracao <= residuo_final else "nao_convergiu"

print(f"=> Resultado: nVolY={nVolY - 2}"
      f" alocacao={end_allocation_time - start_allocation_time:.6f}"
      f" iteracao={end_time - start_iteration_time:.6f}"
      f" numero_iteracao={numero_iteracao}"
      f" residuo={float(residuo_iteracao):.6e}"
      f" estado={estado}")

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
