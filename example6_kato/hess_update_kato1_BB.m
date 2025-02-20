%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = hess_update_kato1_BB(x_new,x_old,y_new,y_old,fun,gj,hess_old)  
% Hessian update for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [2] 
% Experiment 1
% [2] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
%
% Custom Hessian update for experiment 1
% this is an implementarion of the Barzilai and Borwein Hessian update
       
    lambda_min=1e-6;
    lambda_max=1e6;
    [~, fgrad] = fun(x_new);
    [~, fgrad_old] = fun(x_old);
    [~, ggrad] = gj(x_new);
    [~, ggrad_old] = gj(x_old);
    n = length(fgrad);
    
    variation_lagrangian = fgrad - ggrad' * y_new  -fgrad_old + ggrad_old' * y_old;
    dx = x_new-x_old;
    sigma_long = dx'*variation_lagrangian/(dx'*dx);
    %sigma_short = dx'*variation_lagrangian/(variation_lagrangian'*variation_lagrangian);
    sigma=min(lambda_max,max(lambda_min,sigma_long));
    B = sigma * eye(n);

 