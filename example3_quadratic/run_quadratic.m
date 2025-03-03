%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 3: Linear cone contraint Quadratic objective

% Construction of the matrices for the problem
%N=16;
N=20;
%N=35;
%N=50;
vertices=[ones(1,N); zeros(N-1,N)];
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
        matJ((i-1)*N+1:i*N,(j-1)*N+1:j*N)=eye(N);
    end
end
matQ=(N+1)*eye(N^2)-matJ;
% we construct a feasible starting point using the geometry of the problem
x0=10*eye(N);
x0=reshape(x0,N^2,1);
% or we can use the subroutine to look for a feasible point using fdipa
%x0 = searchStartingPoint(N^2,@(x)linear_constraint_fun(x,matA,b),mj);


%my_options = fdipa_options('Display','final','StepTolerance',1e-15,'ParEta',1e-4);
my_options = fdipa_options('Display','final');%,'StepTolerance',1e-15);
[~,fval,~,output1] = fdipaQuad(matQ,cvec,matA,b,x0,mj,[],my_options);
fprintf("Number of Hessian resets: %d \n",output1.hessianresetcount)
fprintf('%d & Exact & %d & %11f & %11.5e & %11f \\\\ %%exact \n',N^2,output1.iterations,fval, output1.firstorderopt, output1.walltime)
 
% % In order to use a different Hessian update we choose the appropriate option
my_options = fdipa_options('Display','final','HessianApproximation','bfgs');%,'StepTolerance',1e-15);
[~,fval,~,output2] = fdipaQuad(matQ,cvec,matA,b,x0,mj,[],my_options);
fprintf("Number of Hessian resets: %d \n",output2.hessianresetcount)
fprintf('%d & BFGS & %d & %11f & %11.5e & %11f \\\\ %%BFGS \n',N^2,output2.iterations,fval, output2.firstorderopt, output2.walltime)

%output for table 10 in the paper
fprintf('%d & %5.2f & %d & %3.2f & %11.5e & %d & %3.2f & %11.5e \\\\  \n',N^2,fval, output2.iterations,output2.walltime, output2.firstorderopt,output1.iterations,output1.walltime, output1.firstorderopt )

clear 'b'  'i' 'j' 'matJ' 'matA' 'matAi' 'matQ' 'cvec' 'mj' ...
    'my_options' 'N' 'vertices' 'x0' 'xmin' 'y0' 'problem' 'x' ...
    'fval' 'exitflag' 'b_update' 'output1' 'output2'

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