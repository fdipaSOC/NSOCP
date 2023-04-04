function [fun,grad_f]=fun_kato1(x,C,d,f)
% Objective function for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
% Experiment 1
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
    x=x(:);
    fun=x'*C*x+f'*x+d'*(x.^4);
    grad_f=2*C*x+f+4*(d.*(x.^3));