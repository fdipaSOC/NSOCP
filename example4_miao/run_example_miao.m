%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 4. Miao-Chen-Ko example of nonlinear convex programs with 
% second-order cone constraints as presented in [2]
% [2] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017

% Example 5.1 of [1]
% $$f(x) = \exp((x_1-3)^2+x_2^2+(x_3-1)^2+(x_4-2)^2+(x_5+1)^2)$$
%
% $$g(x)=x\in\mathcal{K}^5$$
%
disp('Experiment 1, a): using BFGS hessian update');
x0 = [1;0;0;0;0];
my_options = fdipa_options();
[~,fval,~,output] = fdipa(@fun_miao_ex1,x0,@g_miao_ex1,5,[],my_options);

fprintf('Experiment 1 & BFGS & %d & %5.2f & %11.5e & %5.3f \\\\ \\relax %%BFGS \n',output.iterations,fval, output.firstorderopt, output.walltime)


disp('Experiment 1, b): using Modified newton hessian update');
my_options = fdipa_options('HessianApproximation', ...
                    'user-supplied','HessianFcn',@(x_new,x_old,y_new, y_old, fun,gj,hess_old) hess_update_miao_ex1(x_new));
[~,fval,~,output] = fdipa(@fun_miao_ex1,x0,@g_miao_ex1,5,[],my_options);

fprintf('Experiment 1 & Modified Newton & %d & %5.2f & %11.5e & %5.3f \\\\ \\relax %%mod-newton \n',output.iterations,fval, output.firstorderopt, output.walltime);

% Example 5.2 of [1]
% $$f(x) = x_1^2+2x_2^2+2x_1x_2-10x_1-12x_2$$
%
% $$g(x)=\left[\begin{array}{c}
%       8-x_1+3x_2\\
%       3-x_1^2-2x_1+2x_2-x_2^2\\ 
%     \end{array}\right]\in\mathcal{K}^2$$
%
x0 = [1;0];
disp('Experiment 2, a): using BFGS Hessian update');
my_options = fdipa_options();
[~,fval,~,output] = fdipa(@fun_miao_ex2,x0,@g_miao_ex2,[],[],my_options);

fprintf('Experiment 2 & BFGS & %d & %5.2f & %11.5e & %5.3f \\\\ \\relax  %%BFGS \n',output.iterations,fval, output.firstorderopt, output.walltime);

disp('Experiment 2, b): using modified newton Hessian update');
my_options = fdipa_options(...
    'HessianApproximation','user-supplied','HessianFcn',@(x_new,x_old,y_new,y_old,fun,gj,hess_old) hess_update_miao_ex2(x_new,y_new));
[~,fval,~,output] = fdipa(@fun_miao_ex2,x0,@g_miao_ex2,[],[],my_options);
fprintf('Experiment 2 & Modified Newton & %d & %5.2f & %11.5e & %5.3f \\\\ \\relax %%mod-newton \n',output.iterations,fval, output.firstorderopt, output.walltime);

clear 'x0' 'my_options' 'fval' 'output'