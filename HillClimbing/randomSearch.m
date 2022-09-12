% Random Search for optimization
% Jorge Herrada

clear all
close all
clc
plotting = 0; % need the plot? 

% Select function :
% 1 uses (x-2).^2+(y-2).^2
% 2 uses x.*exp(-x.^2 -y.^2)
function_selected = 1; 

if function_selected == 1
    f = @(x,y) (x-2).^2+(y-2).^2; % Function 
    xl = [-5 -5]';  % lower value for x,y coordinates
    xu = [5 5]';    % upper value for x,y coordinates
elseif function_selected == 2
    f = @(x,y) x.*exp(-x.^2 -y.^2);            % Function 
    xl = [-2 -2]';  % lower value for x,y coordinates
    xu = [2 2]';    % upper value for x,y coordinates
end

G = 1000;        % # of iterations

fx_plot = zeros(1,G); % vector to keep record of min value per iteration


fxBest = inf;   % Default value before starting Random Search

xBest = [0 0]';    % Default x,y values before starting Random Search

for g=1:G
    
    xRand = xl+(xu-xl).*rand(2,1); % get random x,y value within range (xl,xu)
    fxRand = f(xRand(1),xRand(2));  % get value of function evaluated in the random coordinate above generated

    % is random better than current best?
    if (fxRand < fxBest)
        % save new best values
        fxBest = fxRand;
        xBest = xRand;
    end

    fx_plot(g) = f(xBest(1),xBest(2)); % save best value on for f(x,y)
%     Plot_Contour(f,[xRand,xBest],xl,xu); % live plotting
end


%  *****************Display results************************

display(['x = ' num2str(xBest(1))]);
display(['y = ' num2str(xBest(2))]);
display(['fx = ' num2str(fxBest)]);


%*******************Plotting*******************************

if plotting
    Plot_Surf(f,xBest,xl,xu);
    figure
    Plot_Contour(f,xBest,xl,xu);
    
    figure
    hold on
    grid on
    plot(fx_plot,'b-','LineWidth',2)
    title("Gráfica de Convergencia")
    xlabel("Iteraciones")
    ylabel("f(x)")
end