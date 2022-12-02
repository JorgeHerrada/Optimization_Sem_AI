% Diferential Evolution using Penalty Function
% DE/best/2/bin

% Jorge Herrada

clear all
close all
clc
plotting = 0; % need the plots? Yes(1), No(0)

% Punto destino 
pd = [-0.6; -0.4];

% define longitud de brazos
a_1 = 0.35;
a_2 = 0.35;
a_3 = 0.25;

% angulos mínimos y máximos de brazos
xl = [-160; -130; -100]*(pi/180); 
xu = [160; 130; 100]*(pi/180);    

% Calculo del punto actual en funcion de los angulos actuales 
p = @(theta_1,theta_2,theta_3) [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1) + a_3*cos(theta_1 + theta_2 + theta_3); ...
                                a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1) + a_3*sin(theta_1 + theta_2 + theta_3)];

% Funcion Objetivo calcula el error entre punto destino y actual
f = @(q) (1/4)*sum((pd-p(q(1),q(2),q(3))).^2); 

% Funcion de penalizacion
fp = @(x,xl,xu) f(x) + 1000*Penalty2(x,xl,xu);

% definicion de parametros iniciales

G = 50;      % # of iterations/generations
D = 3;          % Dimension
N = 10;         % # individuos / particulas

% mejor valor por iteración para grafica de convergencia
f_plot = zeros(1,G);  

% poblacion inicial y su fitness
x = zeros(D,N);         % individuos
fitness = zeros (1,N);  % fitness para cada individuo

% parametros del algoritmo
F = 0.6;    % 0.6, 1.2 - Factor de amplificacion 
CR = 0.9;   % 0.9, 0.6 - Constante de recombinacion


% return
% generamos poblacion inicial
for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);

    fitness(i) = fp(x(:,i),xl,xu);
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
        fu = fp(u,xl,xu);

        % Mutacion es mejor que posicion actual?
        if fu<fitness(i)
            % guardamos mejor posicion y fitness
            x(:,i) = u;
            fitness(i) = fu;
        end

    end

    % guarda indice del mejor valor de la generacion
    [~,xb] = min(fitness);
    f_plot(g) = fitness(xb);
    
    % plot de toda la poblacion de la generacion
    cla
    hold on
    grid on

    for i=1:N
        Dibujar_Manipulador(x(:,i))
    end
    
    plot(pd(1),pd(2),'go','LineWidth',2,'MarkerSize',10)
    drawnow
end




%  *****************Display results************************
x(:,xb)
fitness(xb)
figure
Dibujar_Manipulador(x(:,xb))
plot(pd(1),pd(2),'go','LineWidth',2,'MarkerSize',10)

%*******************Plotting*******************************
if plotting
    
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