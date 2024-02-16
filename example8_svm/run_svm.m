% Requires CVS for matlab [1]
% [1] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined 
% convex programming, version 2.0 beta. 
% http://cvxr.com/cvx, September 2013.
% Download from http://cvxr.com/cvx/download/
% extract and install it by running cvx_startup.m

%Support vector machine example using a Cholesky factorizarion for the
%covariance matrices Sigma_i
% Sigma_i=S_i*S_i';

%The following datasets are obtained from
% Markelle Kelly, Rachel Longjohn, Kolby Nottingham,
% The UCI Machine Learning Repository, https://archive.ics.uci.edu

%load([fileparts(mfilename('fullpath')),'\dataset_Class\breastcancer.mat'])
%load([fileparts(mfilename('fullpath')),'\dataset_Class\diabetes.mat'])
%load([fileparts(mfilename('fullpath')),'\dataset_Class\german_credit.mat'])
load([fileparts(mfilename('fullpath')),'\dataset_Class\splice.mat'])

%data files contains arrays X containing the information asociated with
%each measurement, and an array Y indicating the 
%clasified with the SVM and the vector Y of 1, -1 indicating the
%corresponding clasification of each datapoint.

[m,n]=size(X);

[mu,Mchol_1,Mchol_2]=split_chol(X,Y);
mj=[n+1;n+1;1;1];
%nu=[0.0005;0.0005];
%nu=[0.01;0.005];
nu=[0.005;0.01];
kappa=sqrt(nu./(1-nu));

% We search for a feasible starting point using CVX
cvx_begin quiet
variables w(n) b
subject to
kappa*norm(Mchol_1'*w,2)<=w'*mu(1,:)'+b-1;
kappa*norm(Mchol_2'*w,2)<=-w'*mu(2,:)'-b-1;
cvx_end
x0=[w;b;kappa];

% use BFGS hessian update for the optimzation
myoptions = fdipa_options('HessianApproximation','bfgs');

%experiment 1 : logarithmic objective with weights 1/2, 1/2
disp('Experiment 1: logarithmic objective with weights 1/2, 1/2');
theta=1/2;
fun=@(x)f_svm_CoDo(x,theta);
const=@(x)g_svm(x,mu,Mchol_1,Mchol_2);
[x,fval,exitflag,output]=fdipa(fun,x0,const,mj,[],myoptions);
kappa_opt=x(n+2:end);
Eta_opt=kappa_opt.^2./(1+kappa_opt.^2);
fprintf('fval: %f \n',fval)
fprintf('kappa_opt: %f \n',kappa_opt)
fprintf('Eta_opt: %f \n',Eta_opt)

%Experiment 2: Logatithmic objective with L2 regularization
disp('Experiment 2: Logatithmic objective with L2 regularization and weights C1=C2=1');
%C1=2^(-2); C2=2^(2);
C1=2^(-2); C2=2^(2);
fun=@(x)f_svm_CoDo_L2(x,C1,C2);
const=@(x)g_svm(x,mu,Mchol_1,Mchol_2);
[x,fval,exitflag,output]=fdipa(fun,x0,const,mj,[],myoptions);

%need to understand the results
kappa_opt=x(n+2:end);
Eta_opt=kappa_opt.^2./(1+kappa_opt.^2);
fprintf('fval: %f \n',fval)
fprintf('kappa_opt: %f \n',kappa_opt)
fprintf('Eta_opt: %f \n',Eta_opt)

%what is this for?

% Pred=sign(X*w+b);
% if max(Y)>1
%     find1=find(Y==2);
%     Y(find1)=-1;
% end
% [AUC,Accu,sens,spec]=medi_auc_accu(Pred,Y);
% 