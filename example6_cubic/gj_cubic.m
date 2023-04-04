function [g_fun,grad_g] = gj_cubic(x,A,b)
% Linear constrain for the example of minimization of a cubic function
    x=x(:);
    g_fun = b + A*x;
    grad_g = A;
