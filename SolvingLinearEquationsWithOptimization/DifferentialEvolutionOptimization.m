% Differential Evolution 
% DE/best/2/bin

% Jorge Herrada

clear all
close all
clc


opc = 4; % seleccionar sistema

% definicion del sistemas de ecuaciones
if opc == 1
    A = [2 3; 5 1]; % matriz de escalares
    b = [7;11];     % vector de resultado
    m = numel(b);   % # de ecuaciones
    xl = [-5 -5]';  % lower value for x,y coordinates
    xu = [5 5]';    % upper value for x,y coordinates
elseif opc == 2
    A = [2 3 4; 1 2 3; 5 1 0; 3 4 1]; % matriz de escalares
    b = [19;13;11;13];     % vector de resultado
    m = numel(b);   % # de ecuaciones
    xl = [-5 -5 -5]';  % lower value for x,y coordinates
    xu = [5 5 5]';    % upper value for x,y coordinates
elseif opc == 3
    A = [2 3; 5 4; 2 5; 4 1; 0.5 0.5]; % matriz de escalares
    b = [-5;5;-15;15;0];     % vector de resultado
    m = numel(b);   % # de ecuaciones
    xl = [-5 -5]';  % lower value for x,y coordinates
    xu = [5 5]';    % upper value for x,y coordinates
elseif opc == 4
    A = [3 1 2 2;4 2 2 3;4 3 7 5;1 6 9 8]; % matriz de escalares
    b = [3;5;-3;5];     % vector de resultado
    m = numel(b);       % # de ecuaciones
    xl = [-5 -5 -5 -5]';   % lower value for x,y coordinates
    xu = [5 5 5 5]';      % upper value for x,y coordinates
end

% Function objetivo basada en el error
f = @(x) (1/(2*m))*sum((b - A*x).^2);

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

    % guarda indice del mejor valor de la generacion
%     [~,xb] = min(fitness);
%     f_plot(g) = fitness(xb);
end

[~,xb] = min(fitness);

%  *****************Display results************************
disp(['Mejor valor para x: ']);
x(:,xb)     % mejor valor para x1,x2,...,xn
disp(['Error: ']);
f(x(:,xb))  % error
disp(['Comprobación de Resultado: ']);
A*x(:,xb)   % comprobacion de resultado
