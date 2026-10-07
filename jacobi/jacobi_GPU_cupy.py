import cupy as cp
import sys
import time
# import matplotlib.pyplot as plt

# Conta de um volume interno a partir dos coeficientes e dos quatro vizinhos.
# O CuPy funde as operações da função em um kernel só
@cp.fuse()
def atualiza(aw, pw, ae, pe, a_s, ps, an, pn, bp, ap):
    return (- aw * pw - ae * pe + (- a_s * ps - an * pn) + bp) / ap

# Conta de um volume de face, com três vizinhos, gravada direto no vetor
@cp.fuse()
def atualiza_face(p, a1, p1, a2, p2, a3, p3, bp, ap):
    p[...] = (- a1 * p1 - a2 * p2 - a3 * p3 + bp) / ap

# Início da medição do tempo de alocação
cp.cuda.Device().synchronize()
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
y = cp.linspace(-R - dx / 2, R + dx / 2, nVolY + 2)

nVolX = int(cp.round(L / dx))
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

# Volumes internos: a velocidade de cada linha vale para todas as colunas
ui = u[1:-1, None]

Ap[1:-1, 1:-1] = dx * rho * (- cp.minimum(0, ui) + cp.maximum(0, ui) - min(0, vs) + max(0, vn)) + 4 * Gamma

Aw[1:-1, 1:-1] = -dx * rho * cp.maximum(0, ui) - Gamma
Ae[1:-1, 1:-1] = dx * rho * cp.minimum(0, ui) - Gamma
As[1:-1, 1:-1] = -dx * rho * max(0, vs) - Gamma
An[1:-1, 1:-1] = dx * rho * min(0, vn) - Gamma
Bp[1:-1, 1:-1] = 0

# Solução inicial de phi
phi_new = cp.zeros((nVolY, nVolX))

# Resolução do sistema linear
residuo_iteracao = 1
numero_iteracao = 0
numero_maximo_iteracao = 1000000
residuo_final = 1e-10

# Norma do termo fonte (os cantos fantasmas não fazem parte do sistema)
b = Bp.copy()
b[::nVolY - 1, ::nVolX - 1] = 0
norma_b = cp.linalg.norm(b)

# As fatias são criadas uma vez, antes das iterações. Faces Oeste, Leste, Sul
# e Norte: os volumes a atualizar e os argumentos da conta
faces = [
    (phi_new[1:-1, 0], Ae[1:-1, 0], phi_new[1:-1, 1], As[1:-1, 0], phi_new[:-2, 0], An[1:-1, 0], phi_new[2:, 0], Bp[1:-1, 0], Ap[1:-1, 0]),
    (phi_new[1:-1, -1], Aw[1:-1, -1], phi_new[1:-1, -2], As[1:-1, -1], phi_new[:-2, -1], An[1:-1, -1], phi_new[2:, -1], Bp[1:-1, -1], Ap[1:-1, -1]),
    (phi_new[0, 1:-1], Aw[0, 1:-1], phi_new[0, :-2], Ae[0, 1:-1], phi_new[0, 2:], An[0, 1:-1], phi_new[1, 1:-1], Bp[0, 1:-1], Ap[0, 1:-1]),
    (phi_new[-1, 1:-1], Aw[-1, 1:-1], phi_new[-1, :-2], Ae[-1, 1:-1], phi_new[-1, 2:], As[-1, 1:-1], phi_new[-2, 1:-1], Bp[-1, 1:-1], Ap[-1, 1:-1]),
]

# Volumes internos
volumes = phi_new[1:-1, 1:-1]
argumentos_internos = (Aw[1:-1, 1:-1], phi_new[1:-1, :-2], Ae[1:-1, 1:-1], phi_new[1:-1, 2:], As[1:-1, 1:-1], phi_new[:-2, 1:-1], An[1:-1, 1:-1], phi_new[2:, 1:-1], Bp[1:-1, 1:-1], Ap[1:-1, 1:-1])

print("=> Início das iterações")

# A primeira iteração serve de aquecimento e é descartada
aquecimento = True

cp.cuda.Device().synchronize()
end_allocation_time = time.perf_counter()
start_iteration_time = end_allocation_time

while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
    # Faces Oeste, Leste, Sul e Norte. Os vizinhos que ficam na própria face têm
    # coeficiente nulo, então a conta pode ser gravada direto no vetor
    for argumentos in faces:
        atualiza_face(*argumentos)

    # Volumes internos. A conta inteira é feita com os valores da iteração
    # anterior e só depois gravada no vetor
    volumes[...] = atualiza(*argumentos_internos)

    # Resíduo real do sistema linear: ||b - A*phi|| / ||b||
    res = Bp - Ap * phi_new
    res[:, 1:] -= Aw[:, 1:] * phi_new[:, :-1]
    res[:, :-1] -= Ae[:, :-1] * phi_new[:, 1:]
    res[1:, :] -= As[1:, :] * phi_new[:-1, :]
    res[:-1, :] -= An[:-1, :] * phi_new[1:, :]
    res[::nVolY - 1, ::nVolX - 1] = 0

    residuo_iteracao = cp.linalg.norm(res) / norma_b
    numero_iteracao += 1

    # Aquecimento: a primeira iteração do processo inclui a compilação dos
    # kernels do CuPy. Ela é descartada e o relógio das iterações recomeça
    if aquecimento:
        aquecimento = False
        phi_new[:] = 0
        residuo_iteracao = 1
        numero_iteracao = 0
        cp.cuda.Device().synchronize()
        start_iteration_time = time.perf_counter()

cp.cuda.Device().synchronize()
end_time = time.perf_counter()

estado = "convergiu" if residuo_iteracao <= residuo_final else "nao_convergiu"

print(f"=> Resultado: nVolY={nVolY - 2}"
      f" alocacao={end_allocation_time - start_allocation_time:.6f}"
      f" iteracao={end_time - start_iteration_time:.6f}"
      f" numero_iteracao={numero_iteracao}"
      f" residuo={float(residuo_iteracao):.6e}"
      f" estado={estado}")

# Exibição dos resultados
# x = cp.linspace(0 - dx / 2, L + dx / 2, nVolX)
# X, Y = cp.meshgrid(x, y)
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