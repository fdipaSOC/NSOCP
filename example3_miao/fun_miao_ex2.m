function [fun,grad_f] = fun_miao_ex2(x)
% Objective function for the Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [1]
% Experiment 2
% [1] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017
%
% Function to minimize
% f(x1,x2)=x_1^2+2x_2^2+2x_1*x_2-10*x_1-12*x_2
    x=x(:);
	fun = x(1)^2+2*x(2)^2+2*x(1)*x(2)-10*x(1)-12*x(2);
	grad_f = [2*x(1)+2*x(2)-10;4*x(2)+2*x(1)-12]; 