function [fun,grad_f]=fun_kato2(x,C,d,e,f)
% Objective function for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
% Experiment 2
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
    x=x(:);
    fun=x'*C*x+d'*(x.^4)+e'*(x.^3)+f'*x;
    grad_f=(C + C')*x+4*(d.*(x.^3))+3*(e.*(x.^2))+f;