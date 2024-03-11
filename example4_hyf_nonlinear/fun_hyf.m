function [fun,grad_f] = fun_hyf(x)
% Objective function for Examples as used in [1]
% [1] S. Hayashi, N. Yamashita, and M. Fukushima. A Combined Smoothing
% and Regularization Method for Monotone Second-Order Cone 
% Complementarity Problems. 
% SIAM Journal on Optimization, 15(2):593–615, 2005
    x=x(:);
	fun = exp(x(1)-x(3)) + 3*(2*x(1)-x(2))^4 + sqrt(1+(3*x(2)+5*x(3))^2);
	grad_f = [exp(x(1)-x(3))+24*(2*x(1)-x(2))^3; ...
		-12*(2*x(1)-x(2))^3+3*(3*x(2)+5*x(3))/sqrt(1+(3*x(2)+5*x(3))^2); ...
		-exp(x(1)-x(3))+5*(3*x(2)+5*x(3))/sqrt(1+(3*x(2)+5*x(3))^2)]; 