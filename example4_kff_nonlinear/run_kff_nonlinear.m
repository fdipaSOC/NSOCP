%% Instances of FDIPA for solving nonlinear SOCP problems of the form: 
%  min f(x); 
%  s.t. g^j(x)=A^jx+b^j in K^{mj}, j=1,...,J
%  where
%   f:R^n-->R,  A^j in R^{mj x n}, b^j in R^{mj}   
%
% Examples as used on [1]
%
% Examples are taken of:
% [1] C. Kanzow, I. Ferenczi, and M. Fukushima. On the local convergence 
% of semismooth newton methods for linear and nonlinear second-order 
% cone programs without  strict complementarity. 
% SIAM J. Optim., % 20(1):297-320, 2009.
% https://doi.org/10.1137/060657662
%
% [2] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural 
% network for solving nonlinear convex programs with second-order cone 
% constraints, Information Sciences, Volume 268, 2014, p 255-270, 
% https://doi.org/10.1016/j.ins.2013.10.017
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
my_options = fdipa_options('Display','iter','TolCon',1e-15, ...
    'Maxiter',300,'Hessian','bfgs');

for i = 1:5
   fdipa(@fun_kff,x0{i},[],@g_kff,mj,my_options);
end

x0 = x0{2};

disp('Experiment 1: default options')
fdipa(@fun_kff,x0,[],@g_kff,mj);

disp('Experiment 2: BFGS hessian update')
my_options = fdipa_options('Hessian','bfgs');
fdipa(@fun_kff,x0,[],@g_kff,mj,my_options);

disp('Experiment 3: BFGS hessian update and constraint tolerance 10^(-12)')
my_options = fdipa_options('TolCon',10^(-12),'Hessian','bfgs');
fdipa(@fun_kff,x0,[],@g_kff,mj,my_options);

disp(strcat('Experiment 4: using modified Newton Hessian update and ',...
    'constraint tolerance 10^(-12)'));
my_options = fdipa_options('TolCon',1e-12,'Maxiter',100, ...
    'Hessian','mod-newton','Hessfnc',@hess_update_kff);
fdipa(@fun_kff,x0,[],@g_kff,mj,my_options);

clear 'x0' 'mj' 'my_options' 'i';

      
