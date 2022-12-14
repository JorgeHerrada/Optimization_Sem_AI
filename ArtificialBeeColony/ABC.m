% Artificial Bee Colony ABC

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

% default domain 
xl = [-5 -5]';  % lower value for x,y coordinates
xu = [5 5]';    % upper value for x,y coordinates

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


% parametros del algoritmo
L = 50;     % limite de explotación de fuente de alimento
Pf = 30;    % # de fuentes de alimento / abejas empleadas
Po = N-Pf;  % # de abejas observadoras
lim = zeros(1,Pf);  % Contador de uso de cada fuente de alimento

% mejor valor por iteración para grafica de convergencia
f_plot = zeros(1,G);  

% poblacion inicial y su fitness
x = zeros(D,Pf);            %  fuentes de alimento
fitness = zeros (1,Pf);     % fitness para cada fuente (menor = mejor)
aptitud = zeros(1,Pf);      % aptitud de cada fuente (mayor = mejor)

% generamos poblacion inicial
for i=1:Pf
    % crea solucion random
    x(:,i) = xl+(xu-xl).*rand(D,1);

    % calcula funcion en el punto de la solucion random
    fitness(i) = f(x(1,i),x(2,i));
    
    % calcular aptitud de cada solucion
    if fitness(i)>=0
        aptitud(i) = 1/(1+fitness(i));
    else
        aptitud(i) = 1+abs(fitness(i));
    end

end


for g=1:G

    % ********************** ABEJAS EMPLEADAS *************************
    for i=1:Pf
        % genera k aleatorio distinto de i [1, total_abejas_empleadas]
        k = i;
        while k==i
            k = randi([1 Pf]);
        end
        
        % j toma aleatoriamente una dimension
        j = randi([1 D]);

        % phi toma valor random [-1,1]
        phi = 2*rand()-1;

        % guardar valor de fuente de alimento actual antes de modificarla*
        v = x(:,i);

        % creamos nueva solucion al modificar fuente actual en la dimension 
        % j usando formula del algoritmo
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));

        % evaluamos funcion en nueva solucion y guardamos
        fv = f(v(1),v(2));

        if fv<fitness(i)
            x(:,i) = v;
            fitness(i) = fv;
            lim(i) = 0;
        else
            lim(i) = lim(i) + 1;
        end

        % recalcular aptitud de cada solucion
        if fitness(i)>=0
            aptitud(i) = 1/(1+fitness(i));
        else
            aptitud(i) = 1+abs(fitness(i));
        end
        
    end

    % ********************** ABEJAS OBSERVADORAS *************************
    for i=1:Po

        % seleccionamos indice de abeja trabajadora con base en aptitud
        m = Seleccion(aptitud);

        % genera k aleatorio distinto de i [1, total_abejas_empleadas]
        k = m;
        while k==m
            k = randi([1 Pf]);
        end
        
        % j toma aleatoriamente una dimension
        j = randi([1 D]);

        % phi toma valor random [-1,1]
        phi = 2*rand()-1;

        % guardar valor de fuente de alimento actual antes de modificarla*
        v = x(:,m);

        % creamos nueva solucion al modificar fuente actual en la dimension 
        % j usando formula del algoritmo
        v(j) = x(j,m) + phi*(x(j,m)-x(j,k));

        % evaluamos funcion en nueva solucion y guardamos
        fv = f(v(1),v(2));

        if fv<fitness(m)
            x(:,m) = v;
            fitness(m) = fv;
            lim(m) = 0;
        else
            lim(m) = lim(m) + 1;
        end

        % recalcular aptitud de cada solucion
        if fitness(m)>=0
            aptitud(m) = 1/(1+fitness(m));
        else
            aptitud(m) = 1+abs(fitness(m));
        end
        
    end


    % ********************** ABEJAS EXPLORADORAS *************************
    for i=1:Pf
        if lim(i)>L
            % crea nueva solucion random
            x(:,i) = xl+(xu-xl).*rand(D,1);
        
            % calcula funcion en el punto de la solucion random
            fitness(i) = f(x(1,i),x(2,i));
            
            % calcular aptitud de cada solucion
            if fitness(i)>=0
                aptitud(i) = 1/(1+fitness(i));
            else
                aptitud(i) = 1+abs(fitness(i));
            end

            % reiniciamos contador de explotación
            lim(i) = 0;
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

%% Funciones

function [n] = Seleccion (aptitud)

    aptitud_total = sum(aptitud);
    N = numel(aptitud);

    r = rand();
    P_sum = 0;

    for i=1:N
        P_sum = P_sum + aptitud(i)/aptitud_total;

        if P_sum >= r
            n = i;
            return
        end
    end

    n = N;
end


