clear all
% close all
clc

% define longitud de brazos
a_1 = 0.35;
a_2 = 0.35;
a_3 = 0.25;

% angulo de brazos
% q = [-pi/3; pi/2; pi/3];
q = [-0.3560; 2.2041; -1.0130];

% punto destino
pd = [0.4; 0.4];

% angulo de brazos
theta_1 = q(1);
theta_2 = q(2);
theta_3 = q(3);

% punto actual según el angulo y longitud de brazos
p = [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1) + a_3*cos(theta_1 + theta_2 + theta_3); ...
     a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1) + a_3*sin(theta_1 + theta_2 + theta_3)];

% imprime
cla
Dibujar_Manipulador(q)
plot(pd(1),pd(2),'go','LineWidth',2,'MarkerSize',10)

% punto alcanzado
p
