% Diferential Evolution 
% DE/current to rand/1/exp

% Jorge Herrada

clear all
close all
clc
plotting = 1; % need the plots? Yes(1), No(0)

% 1 uses Griewank
% 2 uses Rastrigin
% 3 uses Drop-Waive
% Select function :
f_selected = 4; 

xl = [-5 -5]';  % lower value for x,y coordinates
xu = [5 5]';    % upper value for x,y coordinates

% Select Function to optimize
if f_selected == 1
    % Griewank Function
    f = @(x,y) ((x.^2/400) + (y.^2/4000)) - (cos(x).*cos(y/sqrt(2))) + 1; 
elseif f_selected == 2
    % Rastrigin Function
    f = @(x,y) 10*2 + (x.^2 - 10*cos(2*pi*x)) + (y.^2 - 10*cos(2*pi*y));
elseif f_selected == 3
    % Drop-Wave function
    f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2)); 
elseif f_selected == 4
    % Mccormick function
    f = @(x,y) sin(x+y)+(x-y).^2-1.5*x+2.5*y+1;
    % define su propio dominio
    xl = [-1.5 -3]';  % lower value for x,y coordinates
    xu = [4 4]';    % upper value for x,y coordinates
end

% definicion de parametros iniciales

G = 150;      % # of iterations/generations
D = 2;          % Dimension
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

    fitness(i) = f(x(1,i),x(2,i));
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

        % Mutamos como DE/current to rand/1/exp
        v = x(:,i) + F*(x(:,r1)-x(:,i)) + F*(x(:,r2)-x(:,r3));

        

        % ************Recombinacion***************
        % recombinacion exp
        u = x(:,i);     % vector de prueba
        j = randi([1,D]);   % randi entre 1,D
        L = 1;      % contador


        while rand()<=CR && L<=D
            u(j) = v(j);
            j = 1 + mod(j,D);
            L = L + 1;
        end

        % ************Seleccion*******************
        % Calculamos fitness del individuo mutado
        fu = f(u(1),u(2));

        % Mutacion es mejor que posicion actual?
        if fu<fitness(i)
            % guardamos mejor posicion y fitness
            x(:,i) = u;
            fitness(i) = fu;
        end

    end

    
    

    % plot de toda la poblacion de la generacion
    Plot_Contour(f,x,xl,xu);
    % guarda indice del mejor valor de la generacion
    [~,xb] = min(fitness);
    f_plot(g) = fitness(xb);
end




%  *****************Display results************************

display(['x = ' num2str(x(1,xb))]);
display(['y = ' num2str(x(2,xb))]);
display(['fx = ' num2str(fitness(xb))]);


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