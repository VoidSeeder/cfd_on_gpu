
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
L     =      2;  % Aproximadamente o comprimento do domínio [m]
D     =   0.01;  % Diâmetro do canal [m]
rho   =    1e3;  % Densidade da água [kg/m3]
mu    = 8.9e-4;  % Viscosidade da água fluido [kg/m3]
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
y       = -R-dx/2:dx:R+dx/2;  % Leva-se em consideração os volumes fantasmas

nVolX   = round(L/dx);        % Adição dos volumes fantasmas
L       = nVolX*dx;           % L deve ser proporcional a D
x       = 0-dx/2:dx:L+dx/2;   % Leva-se em consideração os volumes fantasmas

nVolX   = nVolX+2;
nVolY   = nVolY+2;

[X,Y]   = meshgrid(x,y);    % Matrizes posições
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Campo de velocidade
%--------------------------------------------------------------------------
dPdx = (Pout-Pin)/L;

xw = X-dx/2;
yw = Y;
uw = 1/(4*mu)*(-dPdx)*(R^2-yw.^2);

xe = X+dx/2;
ye = Y;
ue = 1/(4*mu)*(-dPdx)*(R^2-ye.^2);

vs = zeros(size(Y));
vn = zeros(size(Y));

up = 1/(4*mu)*(-dPdx)*(R^2-Y.^2);
vp = zeros(size(Y));
Vp = (up.^2 + vp.^2).^(1/2);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Coeficientes da Eq. da Temperatura
%--------------------------------------------------------------------------
% Face oeste
i       = 1:nVolY;
j       = 1;
Ap(i,j) = 1;
Aw(i,j) = 0;
Ae(i,j) = 1;
As(i,j) = 0;
An(i,j) = 0;
Bp(i,j) = 2*Tin;

% Face leste
i       = 1:nVolY;
j       = nVolX;
Ap(i,j) =  1;
Aw(i,j) = -1;
Ae(i,j) =  0;
As(i,j) =  0;
An(i,j) =  0;
Bp(i,j) =  0;

% Face sul
i       = 1;
j       = 1:nVolX;
Ap(i,j) =  1;
Aw(i,j) =  0;
Ae(i,j) =  0;
As(i,j) =  0;
An(i,j) =  1;
Bp(i,j) =  2*Twall;

% Face norte
i       = nVolY;
j       = 1:nVolX;
Ap(i,j) =  1;
Aw(i,j) =  0;
Ae(i,j) =  0;
As(i,j) =  1;
An(i,j) =  0;
Bp(i,j) =  2*Twall;

% Volumes internos
i       = 2:nVolY-1;
j       = 2:nVolX-1;
Ap(i,j) = dx * rho * ( -min(0,uw(i,j)) +max(0,ue(i,j)) -min(0,vs(i,j)) +max(0,vn(i,j))) + 4*Gamma;
Aw(i,j) = -dx*rho*max(0,uw(i,j)) - Gamma;
Ae(i,j) =  dx*rho*min(0,ue(i,j)) - Gamma;
As(i,j) = -dx*rho*max(0,vs(i,j)) - Gamma;
An(i,j) =  dx*rho*min(0,vn(i,j)) - Gamma;
Bp(i,j) = 0;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Solucao inicial de phi
%--------------------------------------------------------------------------
phi.new     = zeros(size(Y));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Resolucao do sistema linear
%--------------------------------------------------------------------------
residuoIteracao      = 1;
numeroIteracao       = 0;
numeroMaximoIteracao = 10000;
residuoFinal         = 1e-10;

disp('=> Inicio das iteracoes');

while (residuoIteracao > residuoFinal) && (numeroIteracao < numeroMaximoIteracao)
    
    phi.old = phi.new;

    % Face oeste
    for i = 2:nVolY-1
        j = 1;
        phi.new(i,j) = ( - Ae(i,j)*phi.new(i,j+1) - As(i,j)*phi.new(i-1,j) - An(i,j)*phi.new(i+1,j) + Bp(i,j) ) / Ap(i,j);
    end
    
    % Face leste
    for i = 2:nVolY-1
        j = nVolX;
        phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) - As(i,j)*phi.new(i-1,j) - An(i,j)*phi.new(i+1,j) + Bp(i,j) ) / Ap(i,j);
    end
    
    % Face sul
    i   = 1;
    for j = 2:nVolX-1
        phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) - Ae(i,j)*phi.new(i,j+1) - An(i,j)*phi.new(i+1,j) + Bp(i,j) ) / Ap(i,j);
    end
    
    % Face norte
    i = nVolY;
    for j = 2:nVolX-1
        phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) - Ae(i,j)*phi.new(i,j+1) - As(i,j)*phi.new(i-1,j) + Bp(i,j) ) / Ap(i,j);
    end
    
    % Volumes internos
    for i = 2:nVolY-1
        for j = 2:nVolX-1
            phi.new(i,j) = ( -Aw(i,j)*phi.new(i,j-1) - Ae(i,j)*phi.new(i,j+1) - As(i,j)*phi.new(i-1,j) - An(i,j)*phi.new(i+1,j) + Bp(i,j) ) / Ap(i,j);
        end
    end
    
    residuoIteracao = sum(sum(abs(phi.new -phi.old)))/sum(sum(abs(phi.new)));
    numeroIteracao  = numeroIteracao + 1;
    
    disp(['=> Iteracao: ',num2str(numeroIteracao),' Residuo = ',num2str(residuoIteracao)]);
    
%     loglog(numeroIteracao,residuoIteracao,'.k');
%     hold on
%     pause(0.001);
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Exibicao dos resultados
%--------------------------------------------------------------------------
% Campo de velocidade
PLOTVAR.figure	= figure;
hold on
i = 2:nVolY-1; 
j = 2:nVolX-1;
contourf(X(i,j),Y(i,j),Vp(i,j),'EdgeColor','none');
% axis equal; % Ajuste as proporções dos eixos
colorbar('location','EastOutside');
xlim([0 L]);
ylim([-R R]);
title('Campo de velocidade (m/s)');
colorbar('location','EastOutside');
xlabel('x (m)');
ylabel('y (m)');
set(gcf,'color','w');

% Ative shading interp para um gráfico mais suave
shading interp;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Campo de velocidade
PLOTVAR.figure	= figure;
hold on
i = 2:nVolY-1; 
j = 2:nVolX-1;
contourf(X(i,j),Y(i,j),Vp(i,j),'EdgeColor','none');
axis equal; % Ajuste as proporções dos eixos
colorbar('location','EastOutside');
xlim([0 L]);
ylim([-R R]);
title('Campo de velocidade (m/s)');
colorbar('location','EastOutside');
xlabel('x (m)');
ylabel('y (m)');
set(gcf,'color','w');

% Ative shading interp para um gráfico mais suave
shading interp;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Campo de temperatura
i = 2:nVolY-1; 
j = 2:nVolX-1;
PLOTVAR.figure	= figure;
contourf(X(i,j),Y(i,j),phi.new(i,j),'EdgeColor','none');
hold on
% axis equal; % Ajuste as proporções dos eixos
colorbar('location','EastOutside');
xlim([0 L]);
ylim([-R R]);
set(gcf,'color','w');
xlabel('x (m)');
ylabel('y (m)');
title('Campo de Temperatura (ºC)');

colormap(jet);

% Ative shading interp para um gráfico mais suave
shading interp;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Campo de temperatura
i = 2:nVolY-1; 
j = 2:nVolX-1;
PLOTVAR.figure	= figure;
contourf(X(i,j),Y(i,j),phi.new(i,j),'EdgeColor','none');
hold on
axis equal; % Ajuste as proporções dos eixos
colorbar('location','EastOutside');
xlim([0 L]);
ylim([-R R]);
set(gcf,'color','w');
xlabel('x (m)');
ylabel('y (m)');
title('Campo de Temperatura (ºC)');
colormap(jet);

% Ative shading interp para um gráfico mais suave
shading interp;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Perfis de temperatura
%--------------------------------------------------------------------------
PLOTVAR.figure	= figure;
    hold on; % Mantenha o gráfico para o próximo loop

nPos = 5;
% Defina as 5 posições de X para as quais você deseja plotar a variação de temperatura
xPositions = linspace(0+dx,(L-dx),nPos); % substitua x1, x2, x3, x4, x5 pelas suas posições de X

% Cores para as linhas
colors = ['b', 'g', 'r', 'c', 'm']; % azul, verde, vermelho, ciano e magenta

% Para cada posição de X, encontre a temperatura correspondente em todas as posições de Y
for i=1:nPos
    
    % Encontre o índice de X mais próximo da posição atual de X
    distX = x-xPositions(i);
    [~,elemToPlot(i)] = min(abs(distX));
    
    % Pegue a temperatura correspondente em todas as posições de Y
    T2plot(:,i) = phi.new(2:end-1, elemToPlot(i));
        
    % Plote a temperatura versus Y
    plot(T2plot(:,i),y(2:end-1), 'LineStyle', '-', 'Color', colors(i), 'LineWidth', 2);
   
    legendToPlot{i} = [ 'x = ' num2str(xPositions(i),2) ' m']; %#ok<*SAGROW>
end

xlabel('Temperatura (ºC)');
ylabel('y (m)');

legend(legendToPlot{:}); % substitua x1, x2, x3, x4, x5 pelos seus rótulos desejados

grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.5);

%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------');