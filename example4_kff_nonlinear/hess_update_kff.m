function B = hess_update_kff(xkyk) 
% Modified Newton Hessian update for the Kanzow-Ferenczi-Fukushima 
% example of nonlinear convex programs with 
% second-order cone constraints as presented in [1]
% [1] C. Kanzow, I. Ferenczi, and M. Fukushima. On the local convergence 
% of semismooth newton methods for linear and nonlinear second-order 
% cone programs without  strict complementarity. 
% SIAM J. Optim., % 20(1):297-320, 2009.
% https://doi.org/10.1137/060657662
    xkyk = xkyk(:);
    n = 3;
    x= xkyk(1:n);
    B = zeros(n,n);
    
    B(1,1) = exp(x(1)-x(3))+144*(2*x(1)-x(2))^2;
    B(2,1) = -72*(2*x(1)-x(2))^2;
    B(3,1) = -exp(x(1)-x(3));
    
    B(1,2) = -72*(2*x(1)-x(2))^2;
    B(2,2) = 36*(2*x(1)-x(2))^2 + 9*((1+(3*x(2)+5*x(3))^2)^(-1/2) ...
        -(3*x(2)+5*x(3))^2/(1+(3*x(2)+5*x(3))^2)^(3/2));
    B(3,2) = 15*((1+(3*x(2)+5*x(3))^2)^(-1/2) ...
        -(3*x(2)+5*x(3))^2/(1+(3*x(2)+5*x(3))^2)^(3/2));
    
    B(1,3) = B(3,1);
    B(2,3) = B(3,2);
    B(3,3) = exp(x(1)-x(3))+25*((1+(3*x(2)+5*x(3))^2)^(-1/2) ...
        -(3*x(2)+5*x(3))^2/(1+(3*x(2)+5*x(3))^2)^(3/2));
    
    epsilon = min(eig(B));
    if epsilon <= 0.0
    	B = B+(abs(epsilon)+0.1)*eye(n);
    end