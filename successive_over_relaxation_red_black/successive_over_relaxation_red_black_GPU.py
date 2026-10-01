import cupy as cp
import time
# import matplotlib.pyplot as plt

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

nVolY = 400  # Número de volumes na direção y

# Configuração da malha
R = D / 2
dx = D / nVolY
y = cp.linspace(-R - dx / 2, R + dx / 2, nVolY + 2)

nVolX = int(cp.round(L / dx))
L = nVolX * dx
x = cp.linspace(0 - dx / 2, L + dx / 2, nVolX + 2)

nVolX += 2
nVolY += 2

X, Y = cp.meshgrid(x, y)

# Campo de velocidade
dPdx = (Pout - Pin) / L

xw = X - dx / 2
yw = Y
uw = 1 / (4 * mu) * (-dPdx) * (R**2 - yw**2)

xe = X + dx / 2
ye = Y
ue = 1 / (4 * mu) * (-dPdx) * (R**2 - ye**2)

vs = cp.zeros_like(Y)
vn = cp.zeros_like(Y)

up = 1 / (4 * mu) * (-dPdx) * (R**2 - Y**2)
vp = cp.zeros_like(Y)
Vp = cp.sqrt(up**2 + vp**2)

# Coeficientes da equação de temperatura
Ap = cp.zeros((nVolY, nVolX))
Aw = cp.zeros((nVolY, nVolX))
Ae = cp.zeros((nVolY, nVolX))
As = cp.zeros((nVolY, nVolX))
An = cp.zeros((nVolY, nVolX))
Bp = cp.zeros((nVolY, nVolX))

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
Ap[1:-1, 1:-1] = dx * rho * (- cp.minimum(0, uw[1:-1, 1:-1]) + cp.maximum(0, ue[1:-1, 1:-1]) - cp.minimum(0, vs[1:-1, 1:-1]) + cp.maximum(0, vn[1:-1, 1:-1])) + 4 * Gamma

Aw[1:-1, 1:-1] = -dx * rho * cp.maximum(0, uw[1:-1, 1:-1]) - Gamma
Ae[1:-1, 1:-1] = dx * rho * cp.minimum(0, ue[1:-1, 1:-1]) - Gamma
As[1:-1, 1:-1] = -dx * rho * cp.maximum(0, vs[1:-1, 1:-1]) - Gamma
An[1:-1, 1:-1] = dx * rho * cp.minimum(0, vn[1:-1, 1:-1]) - Gamma
Bp[1:-1, 1:-1] = 0

# Solução inicial de phi
phi_new = cp.zeros((nVolY, nVolX))

# Resolução do sistema linear
residuo_iteracao = 1
numero_iteracao = 0
numero_maximo_iteracao = 1000000
residuo_final = 1e-10
omega = 1.0  # Fator de relaxação (1 = Gauss-Seidel red-black)

# Norma do termo fonte (os cantos fantasmas não fazem parte do sistema)
b = Bp.copy()
b[::nVolY - 1, ::nVolX - 1] = 0
norma_b = cp.linalg.norm(b)

print("=> Início das iterações")

startTime = time.time()

while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
    # Volumes vermelhos (i + j par) e depois pretos (i + j ímpar).
    # Os vizinhos de um volume são sempre da outra cor.
    for cor in (0, 1):
        # Face oeste
        i = 2 - cor
        phi_new[i:-1:2, 0] = (1 - omega) * phi_new[i:-1:2, 0] + omega * (
            - Ae[i:-1:2, 0] * phi_new[i:-1:2, 1]
            - As[i:-1:2, 0] * phi_new[i - 1:-2:2, 0]
            - An[i:-1:2, 0] * phi_new[i + 1::2, 0]
            + Bp[i:-1:2, 0]
        ) / Ap[i:-1:2, 0]

        # Face leste
        i = 2 - (cor + nVolX - 1) % 2
        phi_new[i:-1:2, -1] = (1 - omega) * phi_new[i:-1:2, -1] + omega * (
            - Aw[i:-1:2, -1] * phi_new[i:-1:2, -2]
            - As[i:-1:2, -1] * phi_new[i - 1:-2:2, -1]
            - An[i:-1:2, -1] * phi_new[i + 1::2, -1]
            + Bp[i:-1:2, -1]
        ) / Ap[i:-1:2, -1]

        # Face sul
        j = 2 - cor
        phi_new[0, j:-1:2] = (1 - omega) * phi_new[0, j:-1:2] + omega * (
            - Aw[0, j:-1:2] * phi_new[0, j - 1:-2:2]
            - Ae[0, j:-1:2] * phi_new[0, j + 1::2]
            - An[0, j:-1:2] * phi_new[1, j:-1:2]
            + Bp[0, j:-1:2]
        ) / Ap[0, j:-1:2]

        # Face norte
        j = 2 - (cor + nVolY - 1) % 2
        phi_new[-1, j:-1:2] = (1 - omega) * phi_new[-1, j:-1:2] + omega * (
            - Aw[-1, j:-1:2] * phi_new[-1, j - 1:-2:2]
            - Ae[-1, j:-1:2] * phi_new[-1, j + 1::2]
            - As[-1, j:-1:2] * phi_new[-2, j:-1:2]
            + Bp[-1, j:-1:2]
        ) / Ap[-1, j:-1:2]

        # Volumes internos: linhas ímpares e depois linhas pares
        for i in (1, 2):
            j = 2 - (cor + i) % 2
            phi_new[i:-1:2, j:-1:2] = (1 - omega) * phi_new[i:-1:2, j:-1:2] + omega * (
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

    residuo_iteracao = cp.linalg.norm(res) / norma_b
    numero_iteracao += 1
    print(f"=> Iteração: {numero_iteracao}, Resíduo = {residuo_iteracao}")

endTime = time.time()

print(f"Tempo de execução: {endTime - startTime} segundos")

# Exibição dos resultados
# plt.figure()
# plt.contourf(
#     cp.asnumpy(X[1:-1, 1:-1]),
#     cp.asnumpy(Y[1:-1, 1:-1]),
#     cp.asnumpy(phi_new[1:-1, 1:-1]),
#     cmap="jet",
# )
# plt.colorbar(label="Temperatura (°C)")
# plt.title("Campo de Temperatura")
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# plt.axis("equal")
# plt.savefig("campo_temperatura_GPU.png")
# print("=> Gráfico salvo como 'campo_temperatura.png'")
