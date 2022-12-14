% Algoritmo Base

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
y = zeros(D, 1);         % almacena nueva solucion

% parametros del algoritmo
p = 0.8;         % probabilidad
lambda = 1.5;   
sigma2 = (((gamma(1+lambda))/(lambda*gamma((1+lambda)/2)))*((sin((pi*lambda)/2))/(2^((lambda-1)/2))))^(1/lambda);

% generamos poblacion inicial y calculamos fitness
for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);

    fitness(i) = f(x(1,i),x(2,i));
end

% iteramos por cada generacion
for g=1:G

    % seleccionamos al mejor global
    [~,igb] = min(fitness);
    
    % iteramos por cada individuo
    for i=1:N
        
        % polinizacion global o local
        if rand()<p
            % Vuelo de Levy
            u = normrnd(0,sigma2,[D 1]);
            v = normrnd(0,1,[D 1]);
            L = u./(abs(v).^(1/lambda));

            % generacion de nueva solucion con polinizacion global
            y = x(:,i) + L.*(x(:,igb)-x(:,i));
        else
            % definimos j,k aleatorio : i != j != k
            j = i;
            while j==i
                j=randi([1 N]);
            end
            k = j;
            while j==k || j==i
                k=randi([1 N]);
            end

            % generacion de nueva solucion con polinizacion local
            y = x(:,i) + rand()*(x(:,j)-x(:,k));
        end

        % evaluamos nueva solucion generada
        fy = f(y(1),y(2));

        % es mejor que la existente?
        if fy<fitness(i)
            % actualizamos
            x(:,i) = y;
            fitness(i) = fy;
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