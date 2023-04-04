function [g_fun,grad_g] = g_miao_ex1(x)
% Constrain function for the Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [1]
% Experiment 1
% [1] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017
	x=x(:);
	g_fun = x;
	grad_g = eye(5); 