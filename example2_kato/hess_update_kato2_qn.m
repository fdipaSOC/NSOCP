%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB package for nonlinear Second-Order Cone programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hess = hess_update_kato2_qn(x_new,x_old,y_new,y_old,fun,gj,hess_old) 
% Hessian update for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
% Experiment 2
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
%
% Modified Newton update for experiment 2
% F(x)=x'*C*x+sum_{i=1}^n(d_ix_i^4+e_ix_i^3+f_ix_i)
    x = x_new(:);
    y_old = y_old(:);
    nu = x_new - x_old;
    [~,grad_f_old] = fun(x_old);
    [~,grad_f_new] = fun(x_new);
    [~,grad_g_old] = gj(x_old);
    [~,grad_g_new] = gj(x_new);

    wk = grad_f_new - grad_g_new * y_old  - grad_f_old + grad_g_old * y_old;
    if nu'*wk - 0.2*nu' * hess_old * nu >= 0
        theta = 1;
    else
        den = nu'*(hess_old * nu -  wk);
        theta = 0.8 * nu' * hess_old * nu /den;
    end
    uk = theta * wk  + (1-theta) * hess_old * nu;
    hess = hess_old - hess_old * (nu * nu') * hess_old /(nu' * hess_old * nu) + uk * uk'/ (vk' * uk);
 