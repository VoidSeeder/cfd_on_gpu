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

# Vistas em uma dimensão dos mesmos vetores. Nelas, os volumes de uma diagonal
# (i + j constante) ficam igualmente espaçados, de nVolX - 1 em nVolX - 1
# posições, e a diagonal inteira é uma fatia. Os vizinhos são a mesma fatia
# deslocada
Ap_1d = Ap.ravel()
Aw_1d = Aw.ravel()
Ae_1d = Ae.ravel()
As_1d = As.ravel()
An_1d = An.ravel()
Bp_1d = Bp.ravel()
phi_new_1d = phi_new.ravel()
passo = nVolX - 1

# As fatias de todas as diagonais são criadas uma vez, antes das iterações:
# para cada diagonal, os volumes a atualizar, os coeficientes e os vizinhos
diagonais = []

for d in range(2, nVolY + nVolX - 3):
    # Primeira e última linha da diagonal entre os volumes internos
    i_primeiro = max(1, d - (nVolX - 2))
    i_ultimo = min(nVolY - 2, d - 1)

    inicio = d + i_primeiro * passo
    fim = d + i_ultimo * passo + 1

    p = slice(inicio, fim, passo)
    w = slice(inicio - 1, fim - 1, passo)
    e = slice(inicio + 1, fim + 1, passo)
    s = slice(inicio - nVolX, fim - nVolX, passo)
    n = slice(inicio + nVolX, fim + nVolX, passo)

    diagonais.append((phi_new_1d[p], Aw_1d[p], phi_new_1d[w], Ae_1d[p], phi_new_1d[e], As_1d[p], phi_new_1d[s], An_1d[p], phi_new_1d[n], Bp_1d[p], Ap_1d[p]))

print("=> Início das iterações")

end_allocation_time = time.perf_counter()
start_iteration_time = end_allocation_time

while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
    # Um vetor só: cada volume usa os vizinhos já atualizados na iteração
    # atual. A ordem é a do Gauss-Seidel lexicográfico: as faces primeiro e
    # depois os volumes internos. Cada face só depende do volume interno
    # vizinho (os coeficientes entre faces são nulos), então os volumes de
    # uma face não dependem uns dos outros

    # Face oeste
    phi_new[1:-1, 0] = (
        - Ae[1:-1, 0] * phi_new[1:-1, 1]
        - As[1:-1, 0] * phi_new[:-2, 0]
        - An[1:-1, 0] * phi_new[2:, 0]
        + Bp[1:-1, 0]
    ) / Ap[1:-1, 0]

    # Face leste
    phi_new[1:-1, -1] = (
        - Aw[1:-1, -1] * phi_new[1:-1, -2]
        - As[1:-1, -1] * phi_new[:-2, -1]
        - An[1:-1, -1] * phi_new[2:, -1]
        + Bp[1:-1, -1]
    ) / Ap[1:-1, -1]

    # Face sul
    phi_new[0, 1:-1] = (
        - Aw[0, 1:-1] * phi_new[0, :-2]
        - Ae[0, 1:-1] * phi_new[0, 2:]
        - An[0, 1:-1] * phi_new[1, 1:-1]
        + Bp[0, 1:-1]
    ) / Ap[0, 1:-1]

    # Face norte
    phi_new[-1, 1:-1] = (
        - Aw[-1, 1:-1] * phi_new[-1, :-2]
        - Ae[-1, 1:-1] * phi_new[-1, 2:]
        - As[-1, 1:-1] * phi_new[-2, 1:-1]
        + Bp[-1, 1:-1]
    ) / Ap[-1, 1:-1]

    # Volumes internos, por diagonais (i + j = d), em ordem crescente de d.
    # Os vizinhos Oeste e Sul estão na diagonal d - 1, já atualizada na
    # iteração atual, e os vizinhos Leste e Norte na diagonal d + 1, ainda da
    # iteração anterior, como na varredura linha a linha. Os volumes de uma
    # mesma diagonal não dependem uns dos outros
    for p, aw, pw, ae, pe, a_s, ps, an, pn, bp, ap in diagonais:
        p[...] = (- aw * pw - ae * pe - a_s * ps - an * pn + bp) / ap

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
