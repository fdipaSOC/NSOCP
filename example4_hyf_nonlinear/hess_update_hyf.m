function B = hess_update_hyf(x) 
% Modified Newton Hessian update for Examples as used in [1]
% [1] S. Hayashi, N. Yamashita, and M. Fukushima. A Combined Smoothing
% and Regularization Method for Monotone Second-Order Cone 
% Complementarity Problems. 
% SIAM Journal on Optimization, 15(2):593–615, 2005
    n = 3;
    x= x(1:n);
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