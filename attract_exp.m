function [conc,xs,ys] = attract_exp(x)

% computes concentration at point x,y given an exponential decay.  
% source of attractant at [xs,ys]
xs = 5000;
ys = 5000;

conc = exp(-sqrt((xs-x(1))^2+(ys-x(2))^2)/xs);
