function [g_fun,grad_g] = g_kff(x)
% Constrain function for the Kanzow-Ferenczi-Fukushima 
% example of nonlinear convex programs with 
% second-order cone constraints as presented in [1]
% [1] C. Kanzow, I. Ferenczi, and M. Fukushima. On the local convergence 
% of semismooth newton methods for linear and nonlinear second-order 
% cone programs without  strict complementarity. 
% SIAM J. Optim., % 20(1):297-320, 2009.
% https://doi.org/10.1137/060657662
    x=x(:);
	A = [4 6 3;-1 7 -5];
	b = [-1;2];
	g1 = A*x +b;
	g2 = x;
	g_fun = [g1; g2];
	grad_g = [A;eye(3)]; 