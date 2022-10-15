% Genetic Algorithm for optimization
% Jorge Herrada

clear all
% close all
clc
plotting = 1; % need the plot? Yes(1), No(0)

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

Generations = 100;      % # of iterations/generations
Dimension = 2;          % Dimension
N_individuals = 50;     % # of solutions/individuals per iteration
Elite = 10;              % # number of patents in the elite that wont mutate
pMutation = 0.01;        % Mutation Probability

% create a matrix that will contain each individual(as a column vector)
individuals = zeros(Dimension,N_individuals); 

fitness = zeros(1,N_individuals);   % function evaluated on each individual coordinate
aptitude = zeros(1,N_individuals);   % Vector for individual scores

f_plot = zeros(1,Generations);

% ****************Create Initial Population**************
for i=1:N_individuals
    
    % assign random value on x,y within xl,xu range
    individuals(:,i) = xl+(xu-xl).*rand(Dimension,1);
   
end

for g=1:Generations

%     Plot_Contour(f,individuals,xl,xu);

    % ********************Aptitude****************************
    for i=1:N_individuals
        % evaluate function in values generated above
        fitness(i) = f(individuals(1,i),individuals(2,i));
    
        % is solution a positive value? get the aptitude accordingly
        if fitness(i) >= 0
            aptitude(i) = 1/(1+fitness(i));  
        else
            aptitude(i) = 1+abs(fitness(i));
        end
    end

    % save best individual for plotting
    f_plot(g) = min(fitness);
 
    % *********************Selection**************************
    
    % Matrix that will contain children after selection/crossover/mutation
    children = zeros(Dimension,N_individuals - Elite);

    for i=1:2:N_individuals - Elite
        % custom aptitude values just for debugginng purposes
        % aptitude = [10 100 5 10 6 50];
        
        % Get 1st parent index from Roullete
%         p1 = Roulette(aptitude);
%         p1 = Ranking(aptitude);
        p1 = Torneo(aptitude);
        
        % Get 2nd parent index that MUST be different from the 1st one
        p2 = p1; 
        while p1 == p2
%            p2 = Roulette(aptitude); 
%            p2 = Ranking(aptitude); 
            p2 = Torneo(aptitude); 
        end
        
        % get parents coordinates/values for crossover
        parent1 = individuals(:,p1);
        parent2 = individuals(:,p2);
        
        % *********************Crossover**************************
        
        % save parents info in children before crossover
        child1 = parent1;
        child2 = parent2;
        
        % get random crosspoint
        cPoint = randi([1,Dimension]);
        
        % Crossover
        child1(cPoint:Dimension) = parent2(cPoint:Dimension);
        child2(cPoint:Dimension) = parent1(cPoint:Dimension);
        
        children(:,i) = child1;
        children(:,i+1) = child2;
    
    end

    % *********************Mutation**************************

    % Iterate each value(dimension) of each individual and randomly mutate
    for i=1:N_individuals - Elite
        for j=1:Dimension
            
            % pMutation determines if current value of individual mutates
            if rand() < pMutation
                %***1st option for mutation***
%                 children(j,i) = xl(j)+(xu(j)-xl(j))*rand();

                %***2nd option for mutation***
                children(j,i) = children(j,i) + normrnd(0,1);
            end
        end
    end
    
    % ********************Choose Elite**********************

    % sort aptitudes and save the descend indexes
    [~,Indexes] = sort(aptitude,'descend');

    % get value elite indexes (based on best aptitude)
    valueElite = individuals(:,Indexes(1:Elite));

    % individuals die and children become new individuals
    individuals = [children valueElite];

end


% ********************Final aptitude****************************

% Calculate aptitude of all last generation individuals and select the best

    for i=1:N_individuals
        % evaluate function in values generated above
        fitness(i) = f(individuals(1,i),individuals(2,i));
    
        % is solution a positive value? get the aptitude accordingly
        if fitness(i) >= 0
            aptitude(i) = 1/(1+fitness(i));  
        else
            aptitude(i) = 1+abs(fitness(i));
        end
    end

    % Get index of individual with higher aptitude
    [~,i_mejor] = max(aptitude);

    


%  *****************Display results************************

display(['x = ' num2str(individuals(1,i_mejor))]);
display(['y = ' num2str(individuals(2,i_mejor))]);
display(['fx = ' num2str(f(individuals(1,i_mejor),individuals(2,i_mejor)))]);


%*******************Plotting*******************************
if plotting
    Plot_Contour(f,individuals,xl,xu);
    figure
    Plot_Surf(f,individuals,xl,xu);
    
    figure
    hold on
    grid on
    plot(f_plot,'b-','LineWidth',2)
    title("Gráfica de Convergencia")
    xlabel("Iteraciones")
    ylabel("f(x)")
end