function [fun,grad_f] = fun_kff(x)
% Objective function for the Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [1]
% [1] C. Kanzow, I. Ferenczi, and M. Fukushima. On the local convergence 
% of semismooth newton methods for linear and nonlinear second-order 
% cone programs without  strict complementarity. 
% SIAM J. Optim., % 20(1):297-320, 2009.
% https://doi.org/10.1137/060657662
    x=x(:);
	fun = exp(x(1)-x(3)) + 3*(2*x(1)-x(2))^4 + sqrt(1+(3*x(2)+5*x(3))^2);
	grad_f = [exp(x(1)-x(3))+24*(2*x(1)-x(2))^3; ...
		-12*(2*x(1)-x(2))^3+3*(3*x(2)+5*x(3))/sqrt(1+(3*x(2)+5*x(3))^2); ...
		-exp(x(1)-x(3))+5*(3*x(2)+5*x(3))/sqrt(1+(3*x(2)+5*x(3))^2)]; 