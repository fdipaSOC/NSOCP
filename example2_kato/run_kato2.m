%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB package for nonlinear Second-Order Cone programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 2. Kato-Fukushima example of for nonlinear second-order
% cone programs as presented in [1]
% Experiment 2: Non-Linear constrain
% [1] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear 
% second-order cone programs. Optimization Letters 1, 129-144 (2007). 
% https://doi.org/10.1007/s11590-006-0009-2

% set seed
seed = RandStream('mt19937ar','Seed',1);

% choice of cones for the example
%mj= [5;5];
%mj= [5;5;20];
%mj=[5;5;20;20];
%mj=[10;10;20;20;10];
mj=[10;10;20;20;20;20];
mj=[20;20;30;30;20;30;10];
mj=[30;30;40;40;30;30;40];
mj=[40;40;50;50;50;40;30;40];
mj=[50;60;70;70;50;60;60;60];
mj=[60;70;70;70;60;70;60;70;50];
mj=[80;90;90;90;80;100;80;70;60];

n = sum(mj);
nCones=length(mj);
first = ones(nCones,1);
last = mj;
for i=2:nCones
    first(i) = last(i-1)+1;
    last(i) = last(i-1)+mj(i);
end

C=-1 + 2*rand(seed,n,n);
a1=-1 + 2*rand(seed,n,1);
a2=-1 + 2*rand(seed,n,1);
e=-1 + 2*rand(seed,n,1);
f=-1 + 2*rand(seed,n,1);

d=rand(seed,n,1);
b=zeros(n,1);
for i=1:nCones
    b(first(i))= 1;
end
y0 = b;

% % generate random initial data until a feasible starting condition is found
% %default_options = fdipa_options();
% t = 1;
% x0 = (-1+2*rand(seed,n,1))*t;
% %a negative value of lamb_min means that the initial data is unfeasible
% lamb_min = spectral_decomposition(g_kato2_nlin(x0,a1,a2,b),mj);
% 
% while(min(lamb_min) < 0)
%    t = 0.8*t;
%    x0 = (-1+2*rand(seed,n,1))*t;
%    lamb_min = spectral_decomposition(g_kato2_nlin(x0,a1,a2,b),mj);
% end
x0 = zeros(n,1);

%fun_kato2(x0,C,d,e,f)

disp('experiment 2a: Hessian update bfgs');
my_options = fdipa_options('Display', 'final','MaxIterations',10000,'StepTolerance',1e-15);%, ...
    %'StepTolerance',1e-10,'LinearSystemTolerance',1e-4);
[~,fval,~,output]=fdipa(@(x)fun_kato2(x,C,d,e,f),x0,@(x)g_kato2_nlin(x,a1,a2,b), ...
    mj,y0,my_options);
fprintf('[');
fprintf('%g, ', mj(1:end-1));
fprintf('%g]', mj(end));
fprintf(' & %d & %11f & %11.5e & %11f \\\\ \\relax %%BFGS \n',output.iterations,fval, output.firstorderopt, output.cputime)



disp('experiment 2b: Hessian update Mod-Newton');
hess_update = @(x_new, x_old, y_new, y_old, fun, gj, hess_old) hess_update_kato2(x_new,y_new,C,d,e,a1,a2);

my_options = fdipa_options('Display', 'final',...
    'HessianApproximation','used-supplied','HessianFcn', hess_update,...
    'MaxIterations',10000,'StepTolerance',1e-15);%, 'StepTolerance',1e-10,'LinearSystemTolerance',1e-4);   
[~,fval,~,output]=fdipa(@(x)fun_kato2(x,C,d,e,f),x0,@(x)g_kato2_nlin(x,a1,a2,b), ...
    mj,y0,my_options);
fprintf('[');
fprintf('%g, ', mj(1:end-1));
fprintf('%g]', mj(end));
fprintf(' & %d & %11f & %11.5e & %11f \\\\ \\relax %%mod-newton \n',output.iterations,fval, output.firstorderopt, output.cputime)

% disp('experiment 2c: Hessian update Quasi-Newton');
% hess_update = @(x_new, x_old, y_new, y_old, fun, gj, hess_old) hess_update_kato2_qn(x_new,x_old,y_new,y_old,fun,gj,hess_old);
% 
% my_options = fdipa_options('Display', 'final',...
%     'HessianApproximation','used-supplied','HessianFcn', hess_update);    
% [~,fval,~,output]=fdipa(@(x)fun_kato2(x,C,d,e,f),x0,@(x)g_kato2_nlin(x,a1,a2,b), ...
%     mj,y0,my_options);
% fprintf('[');
% fprintf('%g, ', mj(1:end-1));
% fprintf('%g]', mj(end));
% fprintf(' & %d & %11f & %11.5e & %11f \\\\ \\relax %%quasi-newton \n',output.iterations,fval, output.firstorderopt, output.cputime)


clear 'a1' 'a2' 'b' 'C' 'd' 'e' 'f' 'i' 'first' 'lamb_min' 'last' 'mj' ...
    'my_options' 'n' 'nCones' 'seed' 'x0' 'xmin' 'y0' 'x' 'fval' 'exitflag' ...
    'output' 't' 'default_options'

function [fun,grad_f]=fun_kato2(x,C,d,e,f)
% Objective function for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
    x=x(:);
    fun=x'*C*x+d'*(x.^4)+e'*(x.^3)+f'*x;
    grad_f=(C + C')*x+4*(d.*(x.^3))+3*(e.*(x.^2))+f;
end

function [fun_g,grad_g]=g_kato2_nlin(x,a1,a2,b)
% Non-Linear constrain for the Kato-Fukushima example of 
% for nonlinear second-order cone programs as presented in [1] 
    x= x(:);
    n = length(x);
    fun_g=a1.*(exp(x)-1)+b;
    fun_g(n) = fun_g(n)+a2(n)*x(n)*x(1);
    grad_g=zeros(n,n);
    grad_g(n,n) = a1(n)*exp(x(n))+a2(n)*x(1);
    grad_g(1,n) = a2(n)*x(n);
    
    for i=1:(n-1)
        fun_g(i) = fun_g(i) + a2(i)*x(i)*x(i+1);
        grad_g(i,i) = a1(i)*exp(x(i))+a2(i)*x(i+1);
        grad_g(i+1,i) = a2(i)*x(i);
    end
    grad_g = grad_g';
end