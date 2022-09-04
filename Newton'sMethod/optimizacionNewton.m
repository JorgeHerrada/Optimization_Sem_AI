% Newton Method for optimization

clear all
% close all
clc

% Function and it's derivatives
f = @(x)    (20 - 2*x).* (20 - 2*x).* x;    % function
fp = @(x)   12 * x.^2 - 160 * x + 400;      % 1st delivative
fpp = @(x)  24 * x - 160;                   % 2nd derivative

% vector of values for x 
x = 0:0.1:10;

% inital value on "x" to start to approximate
xr = 1;
% # of iterations
N = 10;

% Iterate N times to approximate best value for xr
for i = 1:N
    xr = xr - (fp(xr)/fpp(xr));
end

% Is it a max or a min? 
% Use the 2nd derivative value on xr to find out
if fpp(xr) >= 0
    disp(['min in x = ' num2str(xr) ])
else
    disp(['max in x = ' num2str(xr) ])
    disp(['max volume is = ' num2str(f(xr))])
end

% figure
grid on
hold on
% cla
plot(x,f(x),'b-','LineWidth',2)     % Function plot
plot(x,fp(x),'g--','LineWidth',2)   % 1st derivative plot
plot(x,fpp(x),'r-.','LineWidth',2)  % 2nd derivative plot
plot(xr,f(xr),'r*','LineWidth',2)    % max/min plot

% labels and legends on plot
legend('fx','fp','fpp')
xlabel('x')
ylabel('fx')

