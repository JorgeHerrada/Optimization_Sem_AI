% Evolution Strategies for optimization
% (1 + 1)ES
% Jorge Herrada

clear all
% close all
clc
plotting = 1; % need the plots? Yes(1), No(0)

% Select function :
% 1 uses (x-2).^2+(y-2).^2
% 2 uses x.*exp(-x.^2 -y.^2)
f_selected = 2; 

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
sigma = 0.1;              % Standard Deviation
f_plot = zeros(1,Generations);

% ****************Create Initial Population**************
  
% assign random value on x,y within xl,xu range 
x = xl+(xu-xl).*rand(Dimension,1);
   


for g=1:Generations
    
    % get random vector from the normal distribution with 
    % mean = 0 and standard deviation = sigma
    r = normrnd(0,sigma,[Dimension 1]);

    % generate child by doing parent + random
    y = x + r;

    % is the child generating a lower value in function than parent?
    if f(y(1),y(2)) < f(x(1),x(2))
        % if so, child becomes new patent
        x = y;
    end


    % plot parent and child on each generarion
%     Plot_Contour(f,[x y],xl,xu);
    % saves value of function evaluated in parent values 
    f_plot(g) = f(x(1),x(2));
end







%  *****************Display results************************

display(['x = ' num2str(x(1))]);
display(['y = ' num2str(x(2))]);
display(['fx = ' num2str(f(x(1),x(2)))]);


%*******************Plotting*******************************
if plotting
    Plot_Contour(f,[x y],xl,xu);
    figure
    Plot_Surf(f,[x y],xl,xu);
    
    figure
    hold on
    grid on
    plot(f_plot,'b-','LineWidth',2)
    title("Gráfica de Convergencia")
    xlabel("Iteraciones")
    ylabel("f(x)")
end