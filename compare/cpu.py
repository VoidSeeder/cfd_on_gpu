import numpy as np
import time
# import matplotlib.pyplot as plt

def cpu_script(nVolY):
  start_allocation_time = time.time()
  
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

  # nVolY = 400  # Número de volumes na direção y

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
  numero_maximo_iteracao = 100
  residuo_final = 1e-10

  # print("=> Início das iterações")

  start_iteration_time = time.time()

  while residuo_iteracao > residuo_final and numero_iteracao < numero_maximo_iteracao:
      phi_old = phi_new.copy()

      # # Face oeste
      # j = 0
      # for i in range(1,nVolY-2):
      #     phi_new[i,j] = ( -Ae[i,j]*phi_new[i,j+1] - As[i,j]*phi_new[i-1,j] - An[i,j]*phi_new[i+1,j] + Bp[i,j] ) / Ap[i,j]
              
      # #Face leste
      # j = nVolX - 1
      # for i in range(1,nVolY-2):
      #     phi_new[i,j] = ( -Aw[i,j]*phi_new[i,j-1] - As[i,j]*phi_new[i-1,j] - An[i,j]*phi_new[i+1,j] + Bp[i,j] ) / Ap[i,j]
      
      # # Face sul
      # i = 0
      # for j in range(1,nVolX-2):
      #     phi_new[i,j] = ( -Aw[i,j]*phi_new[i,j-1] - Ae[i,j]*phi_new[i,j+1] - An[i,j]*phi_new[i+1,j] + Bp[i,j] ) / Ap[i,j]
      
      # # Face norte
      # i = nVolY - 1
      # for j in range(1,nVolX-2):
      #     phi_new[i,j] = ( -Aw[i,j]*phi_new[i,j-1] - Ae[i,j]*phi_new[i,j+1] - As[i,j]*phi_new[i-1,j] + Bp[i,j] ) / Ap[i,j]
      
      # Atualiza as bordas Oeste
      phi_new[1:-1, 0] = (
          - Ae[1:-1, 0] * phi_new[1:-1, 1]
          - As[1:-1, 0] * phi_new[0:-2, 0]
          - An[1:-1, 0] * phi_new[2:, 0]
          + Bp[1:-1, 0]
      ) / Ap[1:-1, 0]
      
      # Atualiza as bordas Leste
      phi_new[1:-1, -1] = (
          - Aw[1:-1, -1] * phi_new[1:-1, -2]
          - As[1:-1, -1] * phi_new[0:-2, -1]
          - An[1:-1, -1] * phi_new[2:, -1]
          + Bp[1:-1, -1]
      ) / Ap[1:-1, -1]

      # Atualiza as bordas Sul
      phi_new[0, 1:-1] = (
          - Aw[0, 1:-1] * phi_new[0, 0:-2]
          - Ae[0, 1:-1] * phi_new[0, 2:]
          - An[0, 1:-1] * phi_new[1, 1:-1]
          + Bp[0, 1:-1]
      ) / Ap[0, 1:-1]

      # Atualiza as bordas Norte
      phi_new[-1, 1:-1] = (
          - Aw[-1, 1:-1] * phi_new[-1, 0:-2]
          - Ae[-1, 1:-1] * phi_new[-1, 2:]
          - As[-1, 1:-1] * phi_new[-2, 1:-1]
          + Bp[-1, 1:-1]
      ) / Ap[-1, 1:-1]
      
          
      # Volumes internos
      phi_new[1:-1, 1:-1] = (- Aw[1:-1, 1:-1] * phi_new[1:-1, :-2]
          - Ae[1:-1, 1:-1] * phi_new[1:-1, 2:]
          - As[1:-1, 1:-1] * phi_new[:-2, 1:-1]
          - An[1:-1, 1:-1] * phi_new[2:, 1:-1]
          + Bp[1:-1, 1:-1]) / Ap[1:-1, 1:-1]
              
      residuo_iteracao = np.sum(np.abs(phi_new - phi_old)) / np.sum(np.abs(phi_new))
      numero_iteracao += 1
      # print(f"=> Iteração: {numero_iteracao}, Resíduo = {residuo_iteracao}")

  end_time = time.time()

  # print(f"Tempo de execução: {endTime - startTime} segundos")
  
  return {"allocation": start_iteration_time - start_allocation_time, "iteration": end_time - start_iteration_time}

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
