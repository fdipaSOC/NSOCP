function [fun_g,grad_g]=g_kato1_lin(x,A,b)
% Linear constrain for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
% Experiment 1
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
    x=x(:);
    fun_g=A*x+b;
    grad_g=A;