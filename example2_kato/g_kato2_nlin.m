function [fun_g,grad_g]=g_kato2_nlin(x,a1,a2,b)
% Non-Linear constrain for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
% Experiment 2
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
    x= x(:);
    n = length(x);
    fun_g=a1.*(exp(x)-1)+b;
    fun_g(n) = fun_g(n)+a2(n)*x(n)*x(1);
    grad_g=zeros(n,n);
    grad_g(n,n) = a1(n)*exp(x(n))+a2(n)*x(1);
    grad_g(1,n) = a2(n)*x(n);
    
    for i=1:(n-1)
        fun_g(i) = fun_g(i) + a2(i)*x(i)*x(i+1);
        grad_g(i,i) = a1(i)*exp(x(i))+a2(i)*x(i+1);
        grad_g(i+1,i) = a2(i)*x(i);
    end
    grad_g = grad_g';
