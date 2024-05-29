%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fun,grad_f] = fun_miao_ex1(x)
% Objective function for the Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [2]
% Experiment 1
% [2] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017
%
% Function to minimize
% f(x1,x2,x3,x4,x5)=exp((x_1-3)^2+x_2^2+(x_3-1)^2+(x_4-2)^2+(x_5+1)^2)
    x=x(:);
	fun = exp((x(1)-3)^2+x(2)^2+(x(3)-1)^2+(x(4)-2)^2+(x(5)+1)^2);
	grad_f = 2*fun*(x -[3;0;1;2;-1]); 
