% Diferential Evolution using Penalty Function
% DE/best/2/bin

% Jorge Herrada

clear all
close all
clc
plotting = 0; % need the plots? Yes(1), No(0)

% definicion del sistemas de ecuaciones
opc = 3;

if opc == 1
    A = [2 3; 5 1]; % matriz de escalares
    b = [7;11];     % vector de resultado
    m = numel(b);   % dimension
    xl = [-5 -5]';  % lower value for x,y coordinates
    xu = [5 5]';    % upper value for x,y coordinates
elseif opc == 2
    A = [2 3 4; 1 2 3; 5 1 0; 3 4 1]; % matriz de escalares
    b = [19;13;11;13];     % vector de resultado
    m = numel(b);   % dimension
    xl = [-5 -5 -5]';  % lower value for x,y coordinates
    xu = [5 5 5]';    % upper value for x,y coordinates
elseif opc == 3
    A = [2 3; 5 4; 2 5; 4 1; 0.5 0.5]; % matriz de escalares
    b = [-5;5;-15;15;0];     % vector de resultado
    m = numel(b);   % dimension
    xl = [-5 -5]';  % lower value for x,y coordinates
    xu = [5 5]';    % upper value for x,y coordinates
elseif opc == 4
    A = [2 3; 5 1]; % matriz de escalares
    b = [7;11];     % vector de resultado
    m = numel(b);   % dimension
end

% Function objetivo basada en el error
f = @(x) (1/(2*m))*sum((b - A*x).^2);

% Funciones de penalizacion
% fp = @(x,xl,xu) f(x(1),x(2)) + 1000*Penalty1(x,xl,xu);
% fp = @(x,xl,xu) f(x(1),x(2)) + 1000*Penalty2(x,xl,xu);

% definicion de parametros iniciales

G = 150;      % # of iterations/generations
[~,D] = size(A);          % Dimension
N = 50;         % # individuos / particulas

% mejor valor por iteración para grafica de convergencia
f_plot = zeros(1,G);  

% poblacion inicial y su fitness
x = zeros(D,N);         % individuos
fitness = zeros (1,N);  % fitness para cada individuo

% parametros del algoritmo
F = 0.6;    % 0.6, 1.2 - Factor de amplificacion 
CR = 0.9;   % 0.9, 0.6 - Constante de recombinacion

% generamos poblacion inicial
for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);

    fitness(i) = f(x(:,i));
end


for g=1:G

    for i=1:N
        % **************Mutacion****************

        % obtenemos valores aleatorios para mutaccion, todos distintos
        I = randperm(N);
        I(I == i) = []; % eliminamos i de la permutacion
        % tomamos los primeros elementos de permutacion 
        r1 = I(1);
        r2 = I(2);
        r3 = I(3);
        r4 = I(4);

        % Mutamos como DE/best/2/bin
        [~,best] = min(fitness);
        v = x(:,best) + F*(x(:,r1)-x(:,r2)) + F*(x(:,r3)-x(:,r4));

        

        % ************Recombinacion***************
        u = zeros(D,1);
        k = randi([1,D]);

        for j=1:D 
            if rand()<=CR || j==k
                u(j) = v(j);
            else
                u(j) = x(j,i);
            end
        end

        % ************Seleccion*******************
        % Calculamos fitness del individuo mutado
        fu = f(u);

        % Mutacion es mejor que posicion actual?
        if fu<fitness(i)
            % guardamos mejor posicion y fitness
            x(:,i) = u;
            fitness(i) = fu;
        end

    end

    
    

    % plot de toda la poblacion de la generacion
%     Plot_Contour(f,x,xl,xu);
    % guarda indice del mejor valor de la generacion
    [~,xb] = min(fitness);
    f_plot(g) = fitness(xb);
end




%  *****************Display results************************

x(:,xb)     % mejor valor para x1,x2,...,xn
f(x(:,xb))  % error
A*x(:,xb)   % comprobacion de resultado

%*******************Plotting*******************************
if plotting
    Plot_Contour(f,x(:,xb),xl,xu);
    figure
    Plot_Surf(f,x(:,xb),xl,xu);
    
    figure
    hold on
    grid on
    plot(f_plot,'b-','LineWidth',2)
    title("Gráfica de Convergencia")
    xlabel("Iteraciones")
    ylabel("f(x)")
end

%% Funcion de penalizacion 2 con exponentes
function z = Penalty2(x,xl,xu)
    z = 0;
    D = numel(xl);
    
    for j=1:D
        
        if xl(j)<x(j)
            z = z + 0;
        else
            z = z + (xl(j)-x(j))^2;
        end

        if x(j)<xu(j)
            z = z + 0;
        else
            z = z + (xu(j)-x(j))^2;
        end
    end

end



%% Funcion de penalizacion 1 con suma
function z = Penalty1(x,xl,xu)
    z = 0;
    D = numel(xl);
    
    for j=1:D
        
        if xl(j)<x(j) && x(j)<xu(j)
            z = z + 0;
        else
            z = z + 1;
        end
    end

end




