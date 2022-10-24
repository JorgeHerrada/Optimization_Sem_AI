% Constriction Factor Particle Swarm Optimization (CFPSO)
% Jorge Herrada

clear all
close all
clc
plotting = 1; % need the plots? Yes(1), No(0)

% Select function :
% 1 uses (x-2).^2+(y-2).^2
% 2 uses x.*exp(-x.^2 -y.^2)
f_selected = 3; 

if f_selected == 1
    % Griewank Function
    f = @(x,y) ((x.^2/400) + (y.^2/4000)) - (cos(x).*cos(y/sqrt(2))) + 1; 
elseif f_selected == 2
    % Rastrigin Function
    f = @(x,y) 10*2 + (x.^2 - 10*cos(2*pi*x)) + (y.^2 - 10*cos(2*pi*y));
elseif f_selected == 3
    f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2));
end

xl = [-5 -5]';  % lower value for x,y coordinates
xu = [5 5]';    % upper value for x,y coordinates

% definicion de parametros iniciales

G = 150;      % # of iterations/generations
D = 2;          % Dimension
N = 50;         % # individuos / particulas

% mejor valor por iteración para grafica de convergencia
f_plot = zeros(1,G);  

% ****************Create Initial Population**************

x = zeros(D,N);     % Particulas
v = zeros(D,N);     % velocidad de las particulas
xb = zeros(D,N);    % mejor posición histórica de cada particula

fitness = zeros (1,N);

% parametros del algoritmo

c1 = 2.05;
c2 = 2.05;

phi = c1+c2;
K = 2/(abs(2-phi-sqrt(phi^2-4*phi)));



for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    v(:,i) = 0.3*randn(D,1);
    xb(:,i) = x(:,i);

    fitness(i) = f(x(1,i),x(2,i));

end


for g=1:G

    % Encuentra si la posición actual es mejor que la mejor histórica
    for i=1:N
        % evalua posicion actual
        fx = f(x(1,i),x(2,i));
        % es mejor que la historica?
        if fx < fitness(i)
            % guarda la posición 
            xb(:,i) = x(:,i);
            fitness(i) = fx;
        end

    end

    % guardamos el valor del la mejor particula y su indice 
    [fx_best,I_best] = min(fitness);


    for i=1:N
        % genera nueva velocidad usando inercia(velocidad anterior * w) +  
        % aprendizaje congnitivo (mejor posicion historica vs actual) +
        % aprendizaje social (mejor posición de la población vs actual)
        v(:,i) = K*(v(:,i) + rand()*c1*(xb(:,i)-x(:,i)) + rand()*c2*(xb(:,I_best)-x(:,i)));
        
        % asignamos nueva posicion sumando la velocidad generada
        x(:,i) = x(:,i) + v(:,i);
    end
       

    % plot animation
%     Plot_Contour(f,x,xl,xu);
    % saves value of function evaluated in parent values 
    f_plot(g) = fx_best;
end


%  *****************Display results************************

display(['x = ' num2str(xb(1))]);
display(['y = ' num2str(xb(2))]);
display(['fx = ' num2str(f(xb(1),xb(2)))]);


%*******************Plotting*******************************
if plotting
    Plot_Contour(f,xb,xl,xu);
    figure
    Plot_Surf(f,xb,xl,xu);
    
    figure
    hold on
    grid on
    plot(f_plot,'b-','LineWidth',2)
    title("Gráfica de Convergencia")
    xlabel("Iteraciones")
    ylabel("f(x)")
end