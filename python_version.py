
#==========================================================================
## Exercicio
#==========================================================================
# Autor: Joao Rodrigo Andrade
#
# Tradução linha a linha do código original em MATLAB (octave_version.m).
# Mesmo método (Gauss-Seidel com laços), mesmo critério de parada e mesma
# ordem de varredura. Diferenças inevitáveis da linguagem:
#   - índices começam em 0, e os intervalos 2:n-1 viram range(1, n-1);
#   - phi.new e phi.old viram phi_new e phi_old;
#   - as matrizes de coeficientes precisam ser alocadas antes do uso;
#   - "shading interp" não tem equivalente no contourf do Matplotlib.
#==========================================================================
import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------
plt.close('all')
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
## Parametros de entrada do usuario
#--------------------------------------------------------------------------
L     =      2      # Aproximadamente o comprimento do domínio [m]
D     =   0.01      # Diâmetro do canal [m]
rho   =    1e3      # Densidade da água [kg/m3]
mu    = 8.9e-4      # Viscosidade da água fluido [kg/m3]
Gamma =   0.61/4200 # Coeficiente de condução da água [W/(m K)]
Pin   =      1      # Pressão na entrada [Pa]
Pout  =      0      # Pressão na saída [Pa]
Tin   =     25      # Temperatura de entrada [ºC]
Twall =    100      # Temperatura das paredes [ºC]

nVolY =     15      # Numero de volumes na direção y
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
## Parametros da malha
#--------------------------------------------------------------------------
R       = D/2
dx      = D/nVolY
y       = -R-dx/2 + dx*np.arange(nVolY+2)  # Leva-se em consideração os volumes fantasmas

nVolX   = round(L/dx)                      # Adição dos volumes fantasmas
L       = nVolX*dx                         # L deve ser proporcional a D
x       = 0-dx/2 + dx*np.arange(nVolX+2)   # Leva-se em consideração os volumes fantasmas

nVolX   = nVolX+2
nVolY   = nVolY+2

X, Y    = np.meshgrid(x, y)                # Matrizes posições
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Campo de velocidade
#--------------------------------------------------------------------------
dPdx = (Pout-Pin)/L

xw = X-dx/2
yw = Y
uw = 1/(4*mu)*(-dPdx)*(R**2-yw**2)

xe = X+dx/2
ye = Y
ue = 1/(4*mu)*(-dPdx)*(R**2-ye**2)

vs = np.zeros(Y.shape)
vn = np.zeros(Y.shape)

up = 1/(4*mu)*(-dPdx)*(R**2-Y**2)
vp = np.zeros(Y.shape)
Vp = (up**2 + vp**2)**(1/2)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
## Coeficientes da Eq. da Temperatura
#--------------------------------------------------------------------------
Ap = np.zeros((nVolY, nVolX))
Aw = np.zeros((nVolY, nVolX))
Ae = np.zeros((nVolY, nVolX))
As = np.zeros((nVolY, nVolX))
An = np.zeros((nVolY, nVolX))
Bp = np.zeros((nVolY, nVolX))

# Face oeste
i        = slice(0, nVolY)
j        = 0
Ap[i, j] = 1
Aw[i, j] = 0
Ae[i, j] = 1
As[i, j] = 0
An[i, j] = 0
Bp[i, j] = 2*Tin

# Face leste
i        = slice(0, nVolY)
j        = nVolX-1
Ap[i, j] =  1
Aw[i, j] = -1
Ae[i, j] =  0
As[i, j] =  0
An[i, j] =  0
Bp[i, j] =  0

# Face sul
i        = 0
j        = slice(0, nVolX)
Ap[i, j] =  1
Aw[i, j] =  0
Ae[i, j] =  0
As[i, j] =  0
An[i, j] =  1
Bp[i, j] =  2*Twall

# Face norte
i        = nVolY-1
j        = slice(0, nVolX)
Ap[i, j] =  1
Aw[i, j] =  0
Ae[i, j] =  0
As[i, j] =  1
An[i, j] =  0
Bp[i, j] =  2*Twall

# Volumes internos
i        = slice(1, nVolY-1)
j        = slice(1, nVolX-1)
Ap[i, j] = dx * rho * ( -np.minimum(0, uw[i, j]) +np.maximum(0, ue[i, j]) -np.minimum(0, vs[i, j]) +np.maximum(0, vn[i, j])) + 4*Gamma
Aw[i, j] = -dx*rho*np.maximum(0, uw[i, j]) - Gamma
Ae[i, j] =  dx*rho*np.minimum(0, ue[i, j]) - Gamma
As[i, j] = -dx*rho*np.maximum(0, vs[i, j]) - Gamma
An[i, j] =  dx*rho*np.minimum(0, vn[i, j]) - Gamma
Bp[i, j] = 0
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
## Solucao inicial de phi
#--------------------------------------------------------------------------
phi_new     = np.zeros(Y.shape)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
## Resolucao do sistema linear
#--------------------------------------------------------------------------
residuoIteracao      = 1
numeroIteracao       = 0
numeroMaximoIteracao = 10000
residuoFinal         = 1e-10

print('=> Inicio das iteracoes')

while (residuoIteracao > residuoFinal) and (numeroIteracao < numeroMaximoIteracao):

    phi_old = phi_new.copy()

    # Face oeste
    for i in range(1, nVolY-1):
        j = 0
        phi_new[i, j] = ( - Ae[i, j]*phi_new[i, j+1] - As[i, j]*phi_new[i-1, j] - An[i, j]*phi_new[i+1, j] + Bp[i, j] ) / Ap[i, j]

    # Face leste
    for i in range(1, nVolY-1):
        j = nVolX-1
        phi_new[i, j] = ( -Aw[i, j]*phi_new[i, j-1] - As[i, j]*phi_new[i-1, j] - An[i, j]*phi_new[i+1, j] + Bp[i, j] ) / Ap[i, j]

    # Face sul
    i   = 0
    for j in range(1, nVolX-1):
        phi_new[i, j] = ( -Aw[i, j]*phi_new[i, j-1] - Ae[i, j]*phi_new[i, j+1] - An[i, j]*phi_new[i+1, j] + Bp[i, j] ) / Ap[i, j]

    # Face norte
    i = nVolY-1
    for j in range(1, nVolX-1):
        phi_new[i, j] = ( -Aw[i, j]*phi_new[i, j-1] - Ae[i, j]*phi_new[i, j+1] - As[i, j]*phi_new[i-1, j] + Bp[i, j] ) / Ap[i, j]

    # Volumes internos
    for i in range(1, nVolY-1):
        for j in range(1, nVolX-1):
            phi_new[i, j] = ( -Aw[i, j]*phi_new[i, j-1] - Ae[i, j]*phi_new[i, j+1] - As[i, j]*phi_new[i-1, j] - An[i, j]*phi_new[i+1, j] + Bp[i, j] ) / Ap[i, j]

    residuoIteracao = np.sum(np.sum(np.abs(phi_new -phi_old)))/np.sum(np.sum(np.abs(phi_new)))
    numeroIteracao  = numeroIteracao + 1

    print(f'=> Iteracao: {numeroIteracao} Residuo = {residuoIteracao:.5g}')

#     plt.loglog(numeroIteracao, residuoIteracao, '.k')
#     plt.pause(0.001)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
## Exibicao dos resultados
#--------------------------------------------------------------------------
# Campo de velocidade
plt.figure()
i = slice(1, nVolY-1)
j = slice(1, nVolX-1)
plt.contourf(X[i, j], Y[i, j], Vp[i, j])
# plt.axis('equal') # Ajuste as proporções dos eixos
plt.colorbar(location='right')
plt.xlim([0, L])
plt.ylim([-R, R])
plt.title('Campo de velocidade (m/s)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.gcf().set_facecolor('w')
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Campo de velocidade
plt.figure()
i = slice(1, nVolY-1)
j = slice(1, nVolX-1)
plt.contourf(X[i, j], Y[i, j], Vp[i, j])
plt.axis('equal') # Ajuste as proporções dos eixos
plt.colorbar(location='right')
plt.xlim([0, L])
plt.ylim([-R, R])
plt.title('Campo de velocidade (m/s)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.gcf().set_facecolor('w')
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Campo de temperatura
i = slice(1, nVolY-1)
j = slice(1, nVolX-1)
plt.figure()
plt.contourf(X[i, j], Y[i, j], phi_new[i, j], cmap='jet')
# plt.axis('equal') # Ajuste as proporções dos eixos
plt.colorbar(location='right')
plt.xlim([0, L])
plt.ylim([-R, R])
plt.gcf().set_facecolor('w')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Campo de Temperatura (ºC)')
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Campo de temperatura
i = slice(1, nVolY-1)
j = slice(1, nVolX-1)
plt.figure()
plt.contourf(X[i, j], Y[i, j], phi_new[i, j], cmap='jet')
plt.axis('equal') # Ajuste as proporções dos eixos
plt.colorbar(location='right')
plt.xlim([0, L])
plt.ylim([-R, R])
plt.gcf().set_facecolor('w')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Campo de Temperatura (ºC)')
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Perfis de temperatura
#--------------------------------------------------------------------------
plt.figure()

nPos = 5
# Defina as 5 posições de X para as quais você deseja plotar a variação de temperatura
xPositions = np.linspace(0+dx, (L-dx), nPos)

# Cores para as linhas
colors = ['b', 'g', 'r', 'c', 'm'] # azul, verde, vermelho, ciano e magenta

elemToPlot   = np.zeros(nPos, dtype=int)
T2plot       = np.zeros((nVolY-2, nPos))
legendToPlot = [None]*nPos

# Para cada posição de X, encontre a temperatura correspondente em todas as posições de Y
for i in range(nPos):

    # Encontre o índice de X mais próximo da posição atual de X
    distX = x-xPositions[i]
    elemToPlot[i] = np.argmin(np.abs(distX))

    # Pegue a temperatura correspondente em todas as posições de Y
    T2plot[:, i] = phi_new[1:-1, elemToPlot[i]]

    # Plote a temperatura versus Y
    plt.plot(T2plot[:, i], y[1:-1], linestyle='-', color=colors[i], linewidth=2)

    legendToPlot[i] = f'x = {xPositions[i]:.2g} m'

plt.xlabel('Temperatura (ºC)')
plt.ylabel('y (m)')

plt.legend(legendToPlot)

plt.grid(True, linestyle='--', color='k', alpha=0.5)

#--------------------------------------------------------------------------
print('--------------------------------------------------------------------')

plt.show()
