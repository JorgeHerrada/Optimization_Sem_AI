% Gradient Descent for Optimization
clear all
close all
clc


f = @(x,y) x.* exp(-x.^2 -y.^2);            % Function
Gx = @(x,y) (1-2*x^2)*exp(-x^2-y^2);        % Partial derivative with respect to x
Gy = @(x,y) -2*x*y*exp(-x^2-y^2);           % Partial derivative with respect to y

vector_size = 65; % number of dots within the vectors's range 
iterations = 1000; % number of iterations for algorithm's loop

x_lim = linspace(-2,2,vector_size); % define vector of values for x 
y_lim = linspace(-2,2,vector_size); % define vector of values for y

[X,Y] = meshgrid(x_lim,y_lim); % create matriz with values for X and Y
Z = f(X,Y); % calculate z values based on function f on all x and y coordinates

xy = [-1 -1]';% vector with initial value on x and y

% h minimize the gradient movement (if h < 1) or maximize (if h > 1)
h = 0.01;


% Find min value for 
for i=1:iterations
%     Gx = 2*(xy(1)-2);   % evaluate on x
%     Gy = 2*(xy(2)-2);   % evaluate on y
    G = [Gx(xy(1),xy(2)) Gy(xy(1),xy(2))]';       % gradient

    xy = xy - h * G;    % gets closer to min value with gradient step
end

% Results printed on console
disp(['X = ' num2str(xy(1))])
disp(['Y = ' num2str(xy(2))])
disp(['f(X,Y) = ' num2str(f(xy(1),xy(2)))])



% ******** PLOTING *******
figure
hold on
grid on

surf(X,Y,Z) % plot de la rejilla en 3D
plot3(xy(1),xy(2),f(xy(1),xy(2)),'r*','LineWidth',2,'MarkerSize',10) % plot de un punto cualqueira en 3D
legend({'función','óptimo'},'FontSize',15)

title('Gráfica en 3D','FontSize',15)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
zlabel('f(x,y)','FontSize',15)
view([-20,60]) % estos valores definen la vista 3D del plot

figure
hold on
grid on

contour(X,Y,Z,20) % plot de la rejilla en 2D
plot(xy(1),xy(2),'r*','LineWidth',2,'MarkerSize',10) % plot de un punto cualqueira en 2D
% plot(-0.5,0,'bo','LineWidth',2,'MarkerSize',10) % plot de un punto cualqueira en 2D
legend({'función','óptimo'},'FontSize',15)

title('Gráfica en 2D','FontSize',15)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)