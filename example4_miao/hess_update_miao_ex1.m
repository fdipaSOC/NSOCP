%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = hess_update_miao_ex1(x) 
% Modified Newton Hessian update
% for the Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [2]
% Experiment 1
% [2] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017
    %xkyk = xkyk(:);
    %x= xkyk(1:n);
    n = 5;
    x= x(:);
    B = zeros(n,n);
    f = fun_miao_ex1(x);
    x_aux = [x(1)-3;x(2);x(3)-1;x(4)-2;x(5)+1];
    
    for i=1:n
    	B(i,i) = 2*f*(2*x_aux(i)^2+1);
    end
    for i = 1:(n-1)
    	for j = (i+1):n
    		B(i,j) = 4*f*x_aux(i)*x_aux(j);
    		B(j,i) = 4*f*x_aux(i)*x_aux(j);
    	end
    end
    epsilon = min(eig(B));
    if epsilon <= 0.0
    	B = B+(abs(epsilon)+0.1)*eye(n);
    end
