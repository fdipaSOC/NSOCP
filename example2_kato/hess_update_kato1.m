function B = hess_update_kato1(xkyk,C,d) 
% Hessian update for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
% Experiment 1
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
%
% Custom Hessian update for experiment 1

    xkyk =  xkyk(:);
    n = length(C);
    x= xkyk(1:n);
    B = C+ 12*diag((d.*(x.^2)));
    epsilon = min(eig(B));
    if epsilon <= 0.0
    	B = B+(abs(epsilon)+0.1)*eye(n);
    end

 