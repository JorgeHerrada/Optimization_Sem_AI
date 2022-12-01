% Diferential Evolution using Penalty Function
% DE/best/2/bin

% Jorge Herrada

clear all
close all
clc
plotting = 0; % need the plots? Yes(1), No(0)

% Cargar imagen de qr 
% img_ref = imread("barranca2.jpg");
% img_ref = imread("graduacion1.jpg");
img_ref = imread("moto-rata1.jpg");

% identifica QR y guarda puntos de referencia
[~,~,P] = readBarcode(img_ref,"QR-CODE"); 
[N,M,~] = size(img_ref); % obtener tamaño de imagen

% Carga imagen a desplazar/escalar/rotar
img_des = imread('Jorge.png');
[n,m,~] = size(img_des); % obtener tamaño de imagen


% puntos de referencia del qr
X1 = P(1,:)';
X2 = P(2,:)';
X3 = P(3,:)';

% puntos extremos de imagen a desplazar
x1 = [1 n]';
x2 = [1 1]';
x3 = [m 1]';

xl = [1; 1; -pi; 0];  % lower value for x,y coordinates
xu = [N; M;  pi; 10];    % upper value for x,y coordinates

% Funciones de penalizacion
% fp = @(x,xl,xu) f(x(1),x(2)) + 1000*Penalty1(x,xl,xu);
fp = @(x,xl,xu) f(x(1),x(2)) + 1000*Penalty2(x,xl,xu);

% definicion de parametros iniciales

G = 150;      % # of iterations/generations
D = 4;          % Dimension
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

    % posicion de individuo
    q = x(:,i);
    
    % aplicamos transformacion de similitud a cada punto de la imagen
    % dada su posocion actual
    xp1 = Transformacion_Similitud(q,x1);
    xp2 = Transformacion_Similitud(q,x2);
    xp3 = Transformacion_Similitud(q,x3);
    
    % calculamos distancia entre el error basado en distancia entre 
    % X (puntos de referencia) y xp (puntos transformados)
    e1 = Distancia_Euclidiana(X1,xp1);
    e2 = Distancia_Euclidiana(X2,xp2);
    e3 = Distancia_Euclidiana(X3,xp3);
    
    % definicion de funcion objetivo
    f = (1/6)*(e1^2+e2^2+e3^2);

    fitness(i) = f;
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


        % posicion de individuo
        q = u;
        
        % aplicamos transformacion de similitud a cada punto de la imagen
        % dada su posocion actual
        xp1 = Transformacion_Similitud(q,x1);
        xp2 = Transformacion_Similitud(q,x2);
        xp3 = Transformacion_Similitud(q,x3);
        
        % calculamos distancia entre el error basado en distancia entre 
        % X (puntos de referencia) y xp (puntos transformados)
        e1 = Distancia_Euclidiana(X1,xp1);
        e2 = Distancia_Euclidiana(X2,xp2);
        e3 = Distancia_Euclidiana(X3,xp3);
        
        % definicion de funcion objetivo
        f = (1/6)*(e1^2+e2^2+e3^2);
    
        fu = f + 10000*Penalty2(u,xl,xu);
  
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

fitness(xb)
q = x(:,xb)

Imprimir_Imagenes(q,img_des,img_ref)

return


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

%% Funcion de penalizacion con exponentes
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



%% Funciones para transformacion
function xp = Transformacion_Similitud (qi,xi)
    dx = qi(1);
    dy = qi(2);
    theta = qi(3);
    s = qi(4);
    
    xp = [s*cos(theta) -s*sin(theta); s*sin(theta) s*cos(theta)]*xi + [dx dy]';
end

function e = Distancia_Euclidiana (X,x)
    e = sqrt((X(1)-x(1))^2+(X(2)-x(2))^2);
end

function Imprimir_Imagenes (q,img_des,img_ref)
    dx = q(1);
    dy = q(2);
    theta = q(3);
    s = q(4);
    
    T = [s*cos(theta) -s*sin(theta) dx; s*sin(theta) s*cos(theta) dy; 0 0 1];
    Tp = projective2d(T');
    
    [N,M,~] = size(img_ref);
    [n,m,~] = size(img_des);

    panoramaView = imref2d([N M]);
    Iwarp = imwarp(img_des,Tp,'OutputView',panoramaView);
    Imask = imwarp(true(n,m),Tp,'OutputView',panoramaView);
    
    blender = vision.AlphaBlender('Operation','Binary mask','MaskSource','Input port');
    panorama = step(blender,img_ref,Iwarp,Imask);
    
    imshow(panorama)
end





