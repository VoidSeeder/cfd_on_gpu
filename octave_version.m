%==========================================================================
%% Exercicio
%==========================================================================
% Autor: Joao Rodrigo Andrade
%==========================================================================

%--------------------------------------------------------------------------
clear; clc; close all; format short;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Parametros de entrada do usuario
%--------------------------------------------------------------------------
L     =      2;  % Comprimento do domínio [m]
D     =   0.01;  % Diâmetro do canal [m]
rho   =    1e3;  % Densidade da água [kg/m3]
mu    = 8.9e-4;  % Viscosidade da água [kg/m3]
Gamma =   0.61/4200;  % Coeficiente de condução da água [W/(m K)]
Pin   =      1;  % Pressão na entrada [Pa]   
Pout  =      0;  % Pressão na saída [Pa]
Tin   =     25;  % Temperatura de entrada [ºC]
Twall =    100;  % Temperatura das paredes [ºC]

nVolY =     15;  % Numero de volumes na direção y
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Parametros da malha
%--------------------------------------------------------------------------
R       = D/2;
dx      = D/nVolY;
y       = -R-dx/2:dx:R+dx/2;  % Inclui volumes fantasmas

nVolX   = round(L/dx);        % Inclui volumes fantasmas
L       = nVolX*dx;           % Ajuste de L proporcional a D
x       = 0-dx/2:dx:L+dx/2;   % Inclui volumes fantasmas

nVolX   = nVolX+2;
nVolY   = nVolY+2;

[X,Y]   = meshgrid(x,y);    % Matrizes de posições
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Campo de velocidade
%--------------------------------------------------------------------------
dPdx = (Pout-Pin)/L;

xw = X-dx/2;
yw = Y;
uw = (1/(4*mu))*(-dPdx)*(R^2 - yw.^2);

xe = X+dx/2;
ye = Y;
ue = (1/(4*mu))*(-dPdx)*(R^2 - ye.^2);

vs = zeros(size(Y));
vn = zeros(size(Y));

up = (1/(4*mu))*(-dPdx)*(R^2 - Y.^2);
vp = zeros(size(Y));
Vp = sqrt(up.^2 + vp.^2);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Inicializa coeficientes da equação da temperatura
%--------------------------------------------------------------------------
Ap = zeros(nVolY, nVolX);
Aw = zeros(nVolY, nVolX);
Ae = zeros(nVolY, nVolX);
As = zeros(nVolY, nVolX);
An = zeros(nVolY, nVolX);
Bp = zeros(nVolY, nVolX);

% Face oeste
i = 1:nVolY;
j = 1;
Ap(i,j) = 1;
Ae(i,j) = 1;
Bp(i,j) = 2*Tin;

% Face leste
j = nVolX;
Ap(i,j) = 1;
Aw(i,j) = -1;

% Face sul
i = 1;
j = 1:nVolX;
Ap(i,j) = 1;
An(i,j) = 1;
Bp(i,j) = 2*Twall;

% Face norte
i = nVolY;
Ap(i,j) = 1;
As(i,j) = 1;
Bp(i,j) = 2*Twall;

% Volumes internos
i = 2:nVolY-1;
j = 2:nVolX-1;
Ap(i,j) = dx*rho * ( -min(0,uw(i,j)) + max(0,ue(i,j)) ...
                      -min(0,vs(i,j)) + max(0,vn(i,j)) ) + 4*Gamma;
Aw(i,j) = -dx*rho*max(0,uw(i,j)) - Gamma;
Ae(i,j) =  dx*rho*min(0,ue(i,j)) - Gamma;
As(i,j) = -dx*rho*max(0,vs(i,j)) - Gamma;
An(i,j) =  dx*rho*min(0,vn(i,j)) - Gamma;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Solução inicial
%--------------------------------------------------------------------------
phi.new = zeros(nVolY, nVolX);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Resolução iterativa
%--------------------------------------------------------------------------
residuoIteracao = 1;
numeroIteracao = 0;
numeroMaximoIteracao = 10000;
residuoFinal = 1e-10;

disp('=> Inicio das iteracoes');

while (residuoIteracao > residuoFinal) && (numeroIteracao < numeroMaximoIteracao)
    phi.old = phi.new;

    % Face oeste
    for i = 2:nVolY-1
        j = 1;
        phi.new(i,j) = ( -Ae(i,j)*phi.new(i,j+1) -As(i,j)*phi.new(i-1,j) ...
                          -An(i,j)*phi.new(i+1,j) +Bp(i,j) ) / Ap(i,j);
    end

    % Face leste
    for i = 2:nVolY-1
        j = nVolX;
        phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) -As(i,j)*phi.new(i-1,j) ...
                          -An(i,j)*phi.new(i+1,j) +Bp(i,j) ) / Ap(i,j);
    end

    % Face sul
    i = 1;
    for j = 2:nVolX-1
        phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) -Ae(i,j)*phi.new(i,j+1) ...
                          -An(i,j)*phi.new(i+1,j) +Bp(i,j) ) / Ap(i,j);
    end

    % Face norte
    i = nVolY;
    for j = 2:nVolX-1
        phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) -Ae(i,j)*phi.new(i,j+1) ...
                          -As(i,j)*phi.new(i-1,j) +Bp(i,j) ) / Ap(i,j);
    end

    % Volumes internos
    for i = 2:nVolY-1
        for j = 2:nVolX-1
            phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) -Ae(i,j)*phi.new(i,j+1) ...
                              -As(i,j)*phi.new(i-1,j) -An(i,j)*phi.new(i+1,j) ...
                              +Bp(i,j) ) / Ap(i,j);
        end
    end

    residuoIteracao = sum(sum(abs(phi.new - phi.old)))/sum(sum(abs(phi.new)));
    numeroIteracao  = numeroIteracao + 1;

    disp(['=> Iteracao: ',num2str(numeroIteracao),' Residuo = ',num2str(residuoIteracao)]);
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Plot campo de velocidade
%--------------------------------------------------------------------------
figure;
i = 2:nVolY-1;
j = 2:nVolX-1;
contourf(X(i,j), Y(i,j), Vp(i,j), 'EdgeColor', 'none');
colorbar;
xlim([0 L]);
ylim([-R R]);
title('Campo de velocidade (m/s)');
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'color', 'w');
shading interp;
axis equal;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Plot campo de temperatura
%--------------------------------------------------------------------------
figure;
contourf(X(i,j), Y(i,j), phi.new(i,j), 'EdgeColor', 'none');
colorbar;
xlim([0 L]);
ylim([-R R]);
title('Campo de Temperatura (ºC)');
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'color', 'w');
colormap(jet);
shading interp;
axis equal;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Perfis de temperatura
%--------------------------------------------------------------------------
figure;
hold on;

nPos = 5;
xPositions = linspace(dx, (L - dx), nPos);
colors = ['b', 'g', 'r', 'c', 'm'];

for k = 1:nPos
    [~,elemToPlot] = min(abs(x - xPositions(k)));
    T2plot = phi.new(2:end-1, elemToPlot);
    plot(T2plot, y(2:end-1), 'LineWidth', 2, 'Color', colors(k));
    legendEntries{k} = ['x = ' num2str(xPositions(k),2) ' m'];
end

xlabel('Temperatura (ºC)');
ylabel('y (m)');
legend(legendEntries);
grid on;
title('Perfis de Temperatura em diferentes x');
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.5);

%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------');
