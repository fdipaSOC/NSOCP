function [g_fun,grad_g] = g_hyf(x)
% Constrain function for Examples as used in [1]
% [1] S. Hayashi, N. Yamashita, and M. Fukushima. A Combined Smoothing
% and Regularization Method for Monotone Second-Order Cone 
% Complementarity Problems. 
% SIAM Journal on Optimization, 15(2):593–615, 2005
    x=x(:);
	A = [4 6 3;-1 7 -5];
	b = [-1;2];
	g1 = A*x +b;
	g2 = x;
	g_fun = [g1; g2];
	grad_g = [A;eye(3)]; 