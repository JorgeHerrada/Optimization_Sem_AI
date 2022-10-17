% Evolution Strategies for optimization
% (MU + 1)ES
% Jorge Herrada

clear all
% close all
clc
plotting = 1; % need the plots? Yes(1), No(0)

% Select function :
% 1 uses (x-2).^2+(y-2).^2
% 2 uses x.*exp(-x.^2 -y.^2)
f_selected = 1; 

if f_selected == 1
    f = @(x,y) (x-2).^2+(y-2).^2; % Function 
    xl = [-5 -5]';  % lower value for x,y coordinates
    xu = [5 5]';    % upper value for x,y coordinates
elseif f_selected == 2
    f = @(x,y) x.*exp(-x.^2 -y.^2);            % Function 
    xl = [-2 -2]';  % lower value for x,y coordinates
    xu = [2 2]';    % upper value for x,y coordinates
end

Generations = 1000;      % # of iterations/generations
Dimension = 2;          % Dimension
% sigma = 0.1;            % Standard Deviation
mu = 30;                % # of individuals per generation

% vector that saves min of each iteration and used later to
% plot the convergence of the algorithm
f_plot = zeros(1,Generations); 

% ****************Create Initial Population**************


x = zeros(Dimension,mu+1); 
sigma = zeros(Dimension,mu+1);
fitness = zeros(1,mu+1); 

for i=1:mu
    x(:,i) = xl+(xu-xl).*rand(Dimension,1);
    sigma(:,i) = 0.1.*rand(Dimension,1);
end

% assign random value on x,y within xl,xu range 
% x = xl+(xu-xl).*rand(Dimension,1);
   


for t=1:Generations
    
    r1 = randi([1 mu]);
    r2 = r1;
    
    while r1 == r2
        r2 = randi([1 mu]);
    end

    % Recombinación

    for j=1:Dimension
        if randi([0 1])
            x(j,mu+1) = x(j,r1);
            sigma(j,mu+1) = sigma(j,r1);
        else
            x(j,mu+1) = x(j,r2);
            sigma(j,mu+1) = sigma(j,r2);
        end
    end


    % get random vector from the normal distribution with 
    % mean = 0 and standard deviation = sigma
    r = normrnd(0,sigma(:,mu+1));

    % generate child by doing parent + random
    x(:,mu+1) = x(:,mu+1) + r;
    
    for i=1:mu+1
        fitness(i) = f(x(1,i),x(2,i));
    end

    [~,I] = sort(fitness);

    x = x(:,I);
    sigma = sigma(:,I);
    fitness = fitness(I);


       

    % plot parent and child on each generarion
%     Plot_Contour(f,x,xl,xu);
    % saves value of function evaluated in parent values 
    f_plot(t) = f(x(1,1),x(2,1));
end

xb = x(:,1);

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