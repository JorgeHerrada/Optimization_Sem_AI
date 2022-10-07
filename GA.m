% Genetic Algorithm for optimization
% Jorge Herrada

clear all
close all
clc
plotting = 0; % need the plot? Yes(1), No(0)

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

Generations = 1000;       % # of iterations/generations
Dimension = 2;          % Dimension
N_individuals = 6;         % # of solutions/individuals per iteration

% create a matrix that will contain each individual(as a column vector)
x = zeros(Dimension,N_individuals); 

fitness = zeros(1,N_individuals);   % function evaluated on each individual coordinate
aptitude = zeros(1,N_individuals);   % Vector for individual scores


% ****************Create Initial Population**************
for i=1:N_individuals
    
    % assign random value on x,y within xl,xu range
    x(:,i) = xl+(xu-xl).*rand(Dimension,1);
    
    % evaluate function in values generated above
    fitness(i) = f(x(1,i),x(2,i));

    % is solution a positive value? get the aptitude accordingly
    if fitness(i) >= 0
        aptitude(i) = 1/(1+fitness(i));  
    else
        aptitude(i) = 1+abs(fitness(i));
    end
end

% Just for testing
Plot_Contour(f,x,xl,xu);


% *********************Selection**************************

% custom aptitude values just for debugginng purposes
% aptitude = [10 100 5 10 6 50];

% Get 1st parent index from Roullete
p1 = Roulette(aptitude);

% Get 2nd parent index that MUST be different from the 1st one
p2 = p1; 
while p1 == p2
   p2 = Roulette(aptitude); 
end

% get parents coordinates/values for crossover
parent1 = x(:,p1);
parent2 = x(:,p2);

% *********************Crossover**************************

% save parents info in children before crossover
children1 = parent1;
children2 = parent2;

% get random crosspoint
pc = randi([1,Dimension]);

% Crossover
children1(pc:Dimension) = parent2(pc:Dimension);
children2(pc:Dimension) = parent1(pc:Dimension);


% *********************Mutation**************************







return;









for g=1:Generations
%     Plot_Contour(f,[xBest,yRand],xl,xu);

    

end




%  *****************Display results************************

% display(['x = ' num2str(xBest(1))]);
% display(['y = ' num2str(xBest(2))]);
% display(['fx = ' num2str(fxBest)]);


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