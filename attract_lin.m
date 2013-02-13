function [conc,xs,ys] = attract_lin(x)

% computes concentration at point x,y given an exponential decay.  
% source of attractant at [xs,ys]
xs = [];
ys = 5000;

conc = (1 - abs(ys - x(2))/(ys*2));
