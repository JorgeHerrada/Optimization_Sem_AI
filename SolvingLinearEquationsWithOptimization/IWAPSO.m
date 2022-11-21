% Inertia Weight Particle Swarm Optimization (IWAPSO)
% Jorge Herrada

clear all
% close all
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

G = 150;            % # of iterations/generations
[~,D] = size(A);    % Dimension
N = 50;             % # individuos / particulas

% mejor valor por iteración para grafica de convergencia
f_plot = zeros(1,G);    

% parametros del algoritmo
c1 = 2;
c2 = 2;

w_max = 0.8;
w_min = 0.1;

% ****************Create Initial Population**************

x = zeros(D,N);     % Particulas
v = zeros(D,N);     % velocidad de las particulas
xb = zeros(D,N);    % mejor posición histórica de cada particula

fitness = zeros (1,N);



for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    v(:,i) = 0.3*randn(D,1);
    xb(:,i) = x(:,i);

    fitness(i) = f(x(:,i));

end


for g=1:G

    % Encuentra si la posición actual es mejor que la mejor histórica
    for i=1:N
        % evalua posicion actual
        fx = f(x(:,i));
        % es mejor que la historica?
        if fx < fitness(i)
            % guarda la posición 
            xb(:,i) = x(:,i);
            fitness(i) = fx;
        end

    end

    % guardamos el valor del la mejor particula y su indice 
    [fx_best,I_best] = min(fitness);

    % w se reasigna desde su valor maximo a su minimo, descendiendo en cada
    % iteración
    w = w_max - (g/G)*(w_max-w_min);

    for i=1:N
        % genera nueva velocidad usando inercia(velocidad anterior * w) +  
        % aprendizaje congnitivo (mejor posicion historica vs actual) +
        % aprendizaje social (mejor posición de la población vs actual)
        v(:,i) = w*v(:,i) + rand()*c1*(xb(:,i)-x(:,i)) + rand()*c2*(xb(:,I_best)-x(:,i));
        
        % asignamos nueva posicion sumando la velocidad generada
        x(:,i) = x(:,i) + v(:,i);
    end
       

    % plot animation
%     Plot_Contour(f,x,xl,xu);
    % saves value of function evaluated in parent values 
    f_plot(g) = fx_best;
end

[~,xb] = min(fitness);

%  *****************Display results************************
disp(['Mejor valor para x: ']);
x(:,xb)     % mejor valor para x1,x2,...,xn
disp(['Error: ']);
f(x(:,xb))  % error
disp(['Comprobación de Resultado: ']);
A*x(:,xb)   % comprobacion de resultado