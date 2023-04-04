function B = hess_update_kato2(xkyk,C,d,e,a1,a2) 
% Hessian update for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
% Experiment 2
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129–144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2
%
% Modified Newton update for experiment 2
% F(x)=x'*C*x+sum_{i=1}^n(d_ix_i^4+e_ix_i^3+f_ix_i)
    xkyk = xkyk(:);
    n = length(C);
    x= xkyk(1:n);
    yk = xkyk((n+1):length(xkyk));
    B = C+ 12*diag((d.*(x.^2)))+ 6*diag((e.*x));
    for i=1:n
        B(i,i) = B(i,i) + a1(i)*exp(x(i))*yk(i);
    end
    for i = 1:(n-1)
        for j = (i+1):n
            B(i,j) = B(i,j) + a2(i);
            B(j,i) = B(i,j);
        end
    end
    B(1,n) = B(1,n) + a2(n);
    B(n,1) = B(1,n);
    epsilon = min(eig(B));
    if epsilon <= 0.0
        B = B+(abs(epsilon)+0.1)*eye(n);
    end
 