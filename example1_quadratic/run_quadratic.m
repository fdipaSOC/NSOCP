%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDIPA : feasible directioninterior-point algorithm. See [1] for details
%
% The algorithm follows closely [1] page 1330.
% [1] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) A feasible 
%     direction algorithm for nonlinear second-order cone programs, 
%     Optimization Methods and Software, 34:6, 1322-1341, 
%     DOI: 10.1080/10556788.2018.1506452
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example1: Linear cone contraint Quadratic objective

% Construction of the matrices for the problem
%N=16;
%N=20;
N=35;
%N=50;
vertices=[ones(1,N); zeros(N-1,N)];
x0=10*eye(N);
x0=reshape(x0,N^2,1);
b=-reshape(vertices,N^2,1);
matAi=eye(N);
matA=zeros(N^2);
cvec=zeros(N^2,1);
for i=1:N
    matA((i-1)*N+1:i*N, (i-1)*N+1:i*N)=matAi(:,[2:i,1,(i+1):N]); 
end
mj=N*ones(1,N);
matJ=zeros(N^2);
for i=1:N
    for j=1:N
        matJ((i-1)*N+1:i*N,(j-1)*N+1:j*N)=speye(N);
    end
end
matQ=N*speye(N^2)-matJ;

my_options = fdipa_options('Display','final','StepTolerance',1e-15);
[x,fval,exitflag,output,lambda] = fdipaQuad(matQ,cvec,matA,b,x0,mj,[],my_options);
[~,gradf] = quad_fun(x,matQ,cvec);
[~,gradg] = linear_constraint_fun(x,matA,b);

fprintf('%d & Exact & %d & %11f & %11.5e & %11f \\\\ %%exact \n',N^2,output.iterations,fval, output.firstorderopt, output.cputime)

% % Which is equivalent to an explicit call to fdipa
% 
% problem.objective = @(x) quad_fun(x,matQ,cvec) ;
% problem.x0 = x0;
% problem.y0 = [];
% problem.gj = @(x) linear_constraint_fun(x,matA,b);
% problem.mj = mj;
% my_options.edit('HessianApproximation','user-supplied');
% b_update = @(x_new, x_old, y_new, y_old, fun, gj, hess_old) matQ;
% my_options.edit('HessianFcn',b_update);
% problem.options = my_options;
% [x,fval,exitflag,output] = fdipa(problem);

% In order to use a different Hessian update we choose the appropriate option
my_options = fdipa_options('Display','final','HessianApproximation','bfgs','StepTolerance',1e-15);
[x,fval,exitflag,output] = fdipaQuad(matQ,cvec,matA,b,x0,mj,[],my_options);
fprintf('%d & BFGS & %d & %11f & %11.5e & %11f \\\\ %%BFGS \n',N^2,output.iterations,fval, output.firstorderopt, output.cputime)


clear 'b'  'i' 'j' 'matJ' 'matA' 'matAi' 'matQ' 'cvec' 'mj' ...
    'my_options' 'N' 'vertices' 'x0' 'xmin' 'y0' 'problem' 'x' ...
    'fval' 'exitflag' 'b_update'

function [f,grad_f]=quad_fun(x,Q,c)
    x = x(:);
    f=1/2 *x'*Q*x + c'*x;
    grad_f=Q*x + c;
end

function [g,grad_g]=linear_constraint_fun(x,A,b)
    x= x(:);
    g=A*x+b;
    grad_g=A;
end