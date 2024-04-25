%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB package for nonlinear Second-Order Cone programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 4. Hayashi, Yamashita, Fukushima example of nonlinear 
% convex programs with second-order cone constraints as presented in [1]
% [1] S. Hayashi, N. Yamashita, and M. Fukushima. A Combined Smoothing
% and Regularization Method for Monotone Second-Order Cone 
% Complementarity Problems. 
% SIAM Journal on Optimization, 15(2):593–615, 2005
%
% Non linear problem with different initial conditions
%
% min exp(x(1)-x(3))+3*(2*x(1)-x(2))^4+sqrt(1+(3*x(2)+5*x(3))^2);
% s.t. [4 6 3;-1 7 -5]*x+[-1;2]  in K^2, x in K^3

% Example 4.2 of [1] (also  used in Ex. 5.3 of [2])

x0 = {[1 ;0 ;0] ,  ...
      [1.8860;-0.1890;-0.4081] ,  ...
      [4.3425;0.0875;-0.2332] , ...
      [4.6972;-0.4294;-1.3931] , ...
      [3.2266;-0.7353;-1.5477] , ...
      [3.7282;0.2875;0.2737] };
mj = [2;3]; 

%my_options = fdipa_options('Display','iter','ConstraintTolerance',1e-15);
%for i = 1:5
%   fdipa(@fun_hyf,x0{i},@g_hyf,mj,[],my_options);
%end

x0 = x0{2};

disp('Experiment 1: default options')
[~,fval,~,output] = fdipa(@fun_hyf,x0,@g_hyf,mj);
fprintf('%d & %11f & %11.5e & %11f \\\\ %%BFGS \n',output.iterations,fval, output.firstorderopt, output.cputime)

disp(strcat('Experiment 2: using modified Newton Hessian update'));
hess_update= @(x_new,x_old,y_new,y_old,fun,gj,hess_old) hess_update_hyf(x_new);
my_options = fdipa_options('MaxIterations',100, ...
    'HessianApproximation','user-supplied','HessianFcn',hess_update);
[x,fval,exitflag,output] = fdipa(@fun_hyf,x0,@g_hyf,mj,[],my_options);
fprintf('%d & %11f & %11.5e & %11f \\\\ mod-newtonBFGS \n',output.iterations,fval, output.firstorderopt, output.cputime)

clear 'x0' 'mj' 'my_options' 'i';

      
