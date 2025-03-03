%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = hess_update_kato2(x_new,y_new,C,d,e,a1,a2) 
% Hessian update for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [2] 
% Experiment 2
% [2] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
%
% Modified Newton update for experiment 2
% F(x)=x'*C*x+sum_{i=1}^n(d_ix_i^4+e_ix_i^3+f_ix_i)
    x = x_new(:);
    yk=y_new(:);
    n = length(C);
    B = C+ C' +12*diag((d.*(x.^2)))+ 6*diag((e.*x));
    for i=1:n
        B(i,i) = B(i,i) - a1(i)*exp(x(i))*yk(i);
    end
    for i = 1:(n-1) 
        B(i,i+1) = B(i,i+1) - a2(i)*yk(i);
        B(i+1,i) = B(i,i+1);
    end
    B(1,n) = B(1,n) - a2(n)*yk(n);
    B(n,1) = B(1,n);
     
    try chol(B);
    catch 
        %if cholesky factorization fails indicates that matB is not positive definite.
        B = B+(abs(min(eig(B)))+0.1)*eye(n);
    end
 