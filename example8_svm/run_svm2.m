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

load([fileparts(mfilename('fullpath')),'\dataset_Class\breastcancer.mat'])
%load([fileparts(mfilename('fullpath')),'\dataset_Class\diabetes.mat'])
%load([fileparts(mfilename('fullpath')),'\dataset_Class\german_credit.mat'])
%load([fileparts(mfilename('fullpath')),'\dataset_Class\splice.mat'])

%data files contains arrays X containing the information asociated with
%each measurement, and an array Y indicating the 
%clasified with the SVM and the vector Y of 1, -1 indicating the
%corresponding clasification of each datapoint.

% set seed
%seed = RandStream('mt19937ar','Seed',1);
[m,n]=size(X);

[mu,Mchol_1,Mchol_2]=split_chol(X,Y);
mj=[n+1;n+1;1;1];
%nu=[0.0005;0.0005];
%nu=[0.01;0.005];
%nu=[0.005;0.01];

maxiter = 30;
report = zeros(10,5);

%Experiment 2: Logarithmic objective with L2 regularization
disp('Experiment 2: Logarithmic objective with L2 regularization and weights C1=C2=1, BFGS');
%use BFGS hessian update for the optimzation
myoptions = fdipa_options('NumericalConditioning',1e15);
C1=1; C2=1;
fun=@(x)f_svm_CoDo_L2(x,C1,C2);
const=@(x)g_svm(x,mu,Mchol_1,Mchol_2);

for iter=1:maxiter
     
    nu = 0.1*rand(2,1);
    kappa=sqrt(nu./(1-nu));
    
    % We search for a feasible starting point using CVX
    cvx_begin quiet
    variables w(n) b
    subject to
    kappa*norm(Mchol_1'*w,2)<=w'*mu(1,:)'+b-1;
    kappa*norm(Mchol_2'*w,2)<=-w'*mu(2,:)'-b-1;
    cvx_end
    x0=[w;b;kappa];
    

    [x,fval,~,output]=fdipa(fun,x0,const,mj,[],myoptions);
    kappa_opt=x(n+2:end);
    eta_opt=kappa_opt.^2./(1+kappa_opt.^2);    
    report(iter,:) = [output.iterations, output.cputime, fval,eta_opt(1),eta_opt(2)];
    fprintf('%d & %d & %11f & %11.5e & %11f & %11.5f & %11.5f \\\\ \n',m, output.iterations,fval, output.firstorderopt, output.cputime,eta_opt(1),eta_opt(2))
  
    % kappa_opt=x(n+2:end);
    % Eta_opt=kappa_opt.^2./(1+kappa_opt.^2);
    % fprintf('kappa_opt: %f \n',kappa_opt)
    % fprintf('Eta_opt: %f \n',Eta_opt)
    
    % disp('Experiment 2b: Logarithmic objective with L2 regularization and weights C1=C2=1, mod-newton');
    % b_update = @(x_new,x_old,y_new,y_old,fun,gj,hess_old)  hess_svm2(x_new,y_new,C1,C2,Mchol_1,Mchol_2);
    % myoptions = fdipa_options('HessianApproximation','user-supplied'...
    %     ,'HessianFcn',b_update,'MaxIterations',10000);
    % fun=@(x)f_svm_CoDo_L2(x,C1,C2);
    % const=@(x)g_svm(x,mu,Mchol_1,Mchol_2);
    % [x,fval,~,output]=fdipa(fun,x0,const,mj,[],myoptions);
    % fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',m, output.iterations,fval, output.firstorderopt, output.cputime)
    % kappa_opt=x(n+2:end);
    % Eta_opt=kappa_opt.^2./(1+kappa_opt.^2);
    % fprintf('kappa_opt: %f \n',kappa_opt)
    % fprintf('Eta_opt: %f \n',Eta_opt)
end
max_report = max(report);
min_report = min(report);
mean_report = mean(report);
array2table([max_report; min_report;mean_report] ,...
    'VariableNames',{'iterations','time','fval','eta1', 'eta2'},...
    'RowNames',{'max','min','mean'})

fprintf("{\\bf ds1} & %4.1f & %d & %d & %4.2f & %4.2f & %4.2f & %4.4f & &  & 0 \\\\\n", ...
mean_report(1), min_report(1), max_report(1),mean_report(2), min_report(2), max_report(2),report(1,3))

