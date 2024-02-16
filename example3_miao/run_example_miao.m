% Example 3. Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [1]
% [1] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017

%% Example 5.1 of [1]
% $$f(x) = \exp((x_1-3)^2+x_2^2+(x_3-1)^2+(x_4-2)^2+(x_5+1)^2)$$
%
% $$g(x)=x\in\mathcal{K}^5$$
%
disp('Experiment 1, a): using BFGS hessian update');
x0 = [1;0;0;0;0];
my_options = fdipa_options('ConstraintTolerance',1e-12,'MaxIterations',100,...
    'HessianApproximation','bfgs');
fdipa(@fun_miao_ex1,x0,@g_miao_ex1,5,[],my_options);

disp('Experiment 1, b): using Modified newton hessian update');
my_options = fdipa_options('ConstraintTolerance',1e-12,'MaxIterations',100, ...
    'HessianApproximation','mod-newton','HessianFcn',@hess_update_miao_ex1);
fdipa(@fun_miao_ex1,x0,@g_miao_ex1,5,[],my_options);


%% Example 5.2 of [1]
% $$f(x) = x_1^2+2x_2^2+2x_1x_2-10x_1-12x_2$$
%
% $$g(x)=\left[\begin{array}{c}
%       8-x_1+3x_2\\
%       3-x_1^2-2x_1+2x_2-x_2^2\\ 
%     \end{array}\right]\in\mathcal{K}^2$$
%
x0 = [1;0];
disp('Experiment 2, a): using BFGS Hessian update');
my_options = fdipa_options('ConstraintTolerance',1e-12,'MaxIterations',100,...
    'HessianApproximation','bfgs');
fdipa(@fun_miao_ex2,x0,@g_miao_ex2,[],[],my_options);

disp('Experiment 2, b): using modified newton Hessian update');
my_options = fdipa_options('ConstraintTolerance',1e-12,'MaxIterations',100, ...
    'HessianApproximation','mod-newton','HessianFcn',@hess_update_miao_ex2);
fdipa(@fun_miao_ex2,x0,@g_miao_ex2,[],[],my_options);

clear 'x0' 'my_options' 