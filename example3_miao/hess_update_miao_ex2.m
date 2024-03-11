function B = hess_update_miao_ex2(xk,yk) 
% Modified Newton Hessian update
% for the Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [1]
% Experiment 2
% [1] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017
    %xkyk = xkyk(:);
    %xk=xk(:);
    yk=yk(:);
    n = 2;
    %yk = xkyk((n+1):length(xkyk));
    B = [2 2;2 4] - [2 0;0 2] *yk(2);
    epsilon = min(eig(B));
    if epsilon <= 0.0
    	B = B+(abs(epsilon)+0.1)*eye(n);
    end