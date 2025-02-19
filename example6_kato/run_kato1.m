%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 6. Kato-Fukushima example of for nonlinear second-order
% cone programs as presented in [2]
% Experiment 1: Linear constraint
% [2] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129-144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2

% set seed
seed = RandStream('mt19937ar','Seed',1);

% choice of cones for the example
mj= [5;5];
%mj= [5;5;20];
%mj=[5;5;20;20];
%mj=[10;10;20;20;10];
%mj=[10;10;20;20;20;20];
%mj=[20;20;30;30;20;30;10];
%mj=[30;30;40;40;30;30;40];
%mj=[40;40;50;50;50;40;30;40];
%mj=[50;60;70;70;50;60;60;60];
%mj=[60;70;70;70;60;70;60;70;50];
%mj=[80;90;90;90;80;100;80;70;60];

n = sum(mj);
nCones=length(mj);
first = ones(nCones,1);
last = mj;
for i=2:nCones
    first(i) = last(i-1)+1;
    last(i) = last(i-1)+mj(i);
end

Z=rand(seed,n,n);
C=Z'*Z;
d=rand(seed,n,1);
f=-1 + 2*rand(seed,n,1);
A=2*rand(seed,n,n);
b=zeros(n,1);
for i=1:nCones
    b(first(i))= 1;
end
y0 =b;

x0 = zeros(n,1);

disp('experiment 1a: Hessian update BFGS (default) ');
my_options = fdipa_options('Display', 'final');
[~,fval,~,output] =fdipa(@(x)fun_kato1(x,C,d,f),x0,@(x)g_kato1_lin(x,A,b),mj,y0,my_options);
fprintf("Number of Hessian resets: %d \n",output.hessianresetcount)
%Output for paper
fprintf('~[');
fprintf('%g, ', mj(1:end-1));
fprintf('%g]', mj(end));
fprintf(' & %d & %11f & %11.5e & %11f \\\\ %%BFGS \n',output.iterations,fval, output.firstorderopt, output.walltime)

disp('experiment 1b: Hessian update Modified-Newton ');
hess_update = @(x_new, x_old, y_new, y_old, fun, gj, hess_old) hess_update_kato1(x_new,C,d);
my_options = fdipa_options('Display', 'final',...
    'HessianApproximation','user-supplied','HessianFcn',hess_update);
[~,fval,~,output]=fdipa(@(x)fun_kato1(x,C,d,f),x0,@(x)g_kato1_lin(x,A,b),mj,y0,my_options);
fprintf("Number of Hessian resets: %d \n",output.hessianresetcount)
%Output for paper
fprintf('~[');
fprintf('%g, ', mj(1:end-1));
fprintf('%g]', mj(end));
fprintf(' & %d & %11f & %11.5e & %11f \\\\ %%mod-newton \n',output.iterations,fval, output.firstorderopt, output.walltime)
clear 'A' 'b' 'C' 'd' 'default_options' 'f' 'i' 'first'  'lamb_min' 'last' ...
    'mj' 'my_options' 'n' 'nCones' 'seed' 't' 'x0' 'xmin' 'y0' 'Z'

function [fun,grad_f]=fun_kato1(x,C,d,f)
% Objective function for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [2] 
    x=x(:);
    fun=x'*C*x+f'*x+d'*(x.^4);
    grad_f=2*C*x+f+4*(d.*(x.^3));
end
function [fun_g,grad_g]=g_kato1_lin(x,A,b)
% Linear constrain for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [2] 
    x=x(:);
    fun_g=A*x+b;
    grad_g=A;
end