%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = hess_update_miao_ex2(xk,yk) 
% Modified Newton Hessian update
% for the Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [2]
% Experiment 2
% [2] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017
    %xkyk = xkyk(:);
    %xk=xk(:);
    yk=yk(:);
    n = 2;
    %yk = xkyk((n+1):length(xkyk));
    B = [2 2;2 4] - [2 0;0 2] *yk(2);
    
    try chol(B);
    catch 
        %if cholesky factorization fails indicates that matB is not positive definite.
        B = B+(abs(min(eig(B)))+0.1)*eye(n);
    end