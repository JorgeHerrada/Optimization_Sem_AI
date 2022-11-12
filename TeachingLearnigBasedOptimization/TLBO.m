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

G = 100;      % # of iterations/generations
D = 2;          % Dimension
N = 30;         % # individuos / particulas

% mejor valor por iteración para grafica de convergencia
f_plot = zeros(1,G);  

% poblacion inicial y su fitness
x = zeros(D,N);         % individuos
fitness = zeros (1,N);  % fitness para cada individuo

% parametros del algoritmo


% generamos poblacion inicial y evaluamos en funcion para obtener fitness
for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(1,i),x(2,i));

end

% Itera Generaciones
for g=1:G

    % Itera individuos
    for i=1:N
        % ******************* Enseñanza ************************
        % determinar mejor actual como maestro 
        [~,t] = min(fitness);
        Tf = randi([1,2]);      % factor de enseñanza
        c = zeros(D,1);         % nueva solucion

        % Itera Dimension (materia)
        for j=1:D
            % encuentra la media en la dimension J de los individuos x
            x_mean = mean(x(j,:));

            % Genera nueva solucion para dimension j
            c(j) = x(j,i)+rand()*(x(j,t)-(Tf*x_mean));
        end

        % evaluamos nueva solucion en funcion objetivo
        fc = f(c(1),c(2));

        % la nueva solucion es mejor?
        if fc<fitness(i)
            % si es asi, guardamos 
            x(:,i) = c;
            fitness(i) = fc;
        end


        % ******************* Aprendizaje ************************ 
        % generamos aleatorio != i [1,N]
        k = i;
        while k == i
            k = randi([1,N]);
        end

        % es mejor la solucion actual o la del otro compañero
        if fitness(i) < fitness(k)
            % si es mejor la actual nos alejamos del compañero
            for j=1:D
                c(j) = x(j,i)+rand()*(x(j,i)-x(j,k));
            end
        else
            % si es mejor el nuevo, nos acercamos al compañero
            for j=1:D
                c(j) = x(j,i)+rand()*(x(j,k)-x(j,i));
            end
        end

        % evaluamos nueva solucion en funcion objetivo
        fc = f(c(1),c(2));

        % la nueva solucion es mejor?
        if fc<fitness(i)
            % si es asi, guardamos 
            x(:,i) = c;
            fitness(i) = fc;
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