% Random Mutation Hill Climbing for optimization
% Jorge Herrada

clear all
close all
clc
plotting = 0; % need the plot? 

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

G = 1000;        % # of iterations
D = 2;          % Dimension

fx_plot = zeros(1,G); % vector to keep record of min value per iteration


xBest = xl+(xu-xl).*rand(D,1); % get random x,y value within range (xl,xu)

for g=1:G
    fxBest = f(xBest(1),xBest(2));

    yRand = xBest; % xBest backup to get randomly mutated 

    % Choose random axis to move on 
    j = randi([1,D]);  

    r = rand(); % random (0,1)

    % get new random value for j's axis of yRand
    yRand(j) = xl(j)+(xu(j)-xl(j))*r; 

    % evaluate function with an axis modified (mutated)
    fyRand = f(yRand(1),yRand(2));

    % is the new random mutation better than current best?
    if fyRand < fxBest 
        xBest = yRand; % new best
    end

    fx_plot(g) = f(xBest(1),xBest(2)); % save best value on for f(x,y)
%     Plot_Contour(f,[xBest,yRand],xl,xu);
end

%  *****************Display results************************

display(['x = ' num2str(xBest(1))]);
display(['y = ' num2str(xBest(2))]);
display(['fx = ' num2str(fxBest)]);


%*******************Plotting*******************************
if plotting
    Plot_Contour(f,xBest,xl,xu);
    figure
    Plot_Surf(f,xBest,xl,xu);
    
    figure
    hold on
    grid on
    plot(fx_plot,'b-','LineWidth',2)
    title("Gráfica de Convergencia")
    xlabel("Iteraciones")
    ylabel("f(x)")
end