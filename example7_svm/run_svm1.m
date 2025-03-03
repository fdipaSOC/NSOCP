%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%An application to Robust Support Vector Machines

% Requires CVS for matlab [2]
% [2] Michael Grant and Stephen Boyd. CVX: Matlab software for disciplined 
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
load([fileparts(mfilename('fullpath')),'\dataset_Class\breastcancer.mat']); label = 'ds1';
%load([fileparts(mfilename('fullpath')),'\dataset_Class\diabetes.mat']); label = 'ds2';
%load([fileparts(mfilename('fullpath')),'\dataset_Class\german_credit.mat']); label = 'ds3';
%load([fileparts(mfilename('fullpath')),'\dataset_Class\splice.mat']); label = 'ds4';

%data files contains arrays X containing the information asociated with
%each measurement, and an array Y indicating the 
%clasified with the SVM and the vector Y of 1, -1 indicating the
%corresponding clasification of each datapoint.

% set seed
seed = RandStream('mt19937ar','Seed',1);
[m,n]=size(X);

[mu,Mchol_1,Mchol_2]=split_chol(X,Y);
mj=[n+1;n+1;1;1];
theta=1/2;
maxiter = 30;
report = zeros(maxiter,8);
Prediction_X = zeros(m,maxiter);
AUC = zeros(maxiter,1);
Accu = zeros(maxiter,1);

%%experiment 1 : logarithmic objective with weights 1/2, 1/2
disp('Experiment 1: logarithmic objective with weights 1/2, 1/2, BFGS');
% use BFGS hessian update for the optimzation
myoptions = fdipa_options('MaxIterations',10000,'Display','final');

fun=@(x)f_svm_CoDo(x,theta);
const=@(x)g_svm(x,mu,Mchol_1,Mchol_2);

for iter=1:maxiter
    nu = 0.1*rand(seed,2,1);
    kappa=sqrt(nu./(1-nu));

    %look for feasible starting point using an auxiliary problem
    %in order to find different feasible points, take a different 
    % initial guess where kappa is chosen randomly
    xguess = [zeros(n+1,1);kappa];
    [x0,~] = searchStartingPoint(n+3,const,mj,xguess);
   
    [x,fval,~,output]=fdipa(fun,x0,const,mj,[],myoptions);
    kappa_opt=x(n+2:end);
    eta_opt=kappa_opt.^2./(1+kappa_opt.^2);    
    report(iter,:) = [output.iterations, output.walltime, fval,output.firstorderopt , output.compslack,output.constrviolation,eta_opt(1),eta_opt(2)];
    fprintf('%d & %d & %11f & %11.5e & %11f & %11.5f & %11.5f \\\\ \n',m, output.iterations,fval, output.firstorderopt, output.walltime,eta_opt(1),eta_opt(2))
    Prediction_X(:,iter)=sign(X*x(1:n)+x(n+1));
    [AUC(iter),Accu(iter)]=medi_auc_accu(Prediction_X(:,iter),Y);
end

mean_AUC=mean(AUC);
max_Accu=max(Accu);
min_Accu=min(Accu);
mean_Accu=mean(Accu);

max_report = max(report);
min_report = min(report);
mean_report = mean(report);
summary=array2table([max_report max_Accu ; min_report min_Accu;mean_report  mean_Accu] ,...
    'VariableNames',{'iterations','time','fval','norm_lag','comp_slack','feasibility','eta1', 'eta2','Accu'},...
    'RowNames',{'max','min','mean'});
disp(summary);

fprintf("{\\bf " + label + "} & %4.1f & %d & %d & %4.2f & %4.2f & %4.2f & %4.4f & %4.4f & %11.5e & %11.5e & %4.4f\\\\\n", ...
mean_report(1), min_report(1), max_report(1),mean_report(2), min_report(2), max_report(2),mean_report(1,7),mean_report(1,8),mean_report(1,4),mean_report(1,5),mean_Accu);

clear 'Accu' 'AUC' 'const' 'eta_opt' 'fun' 'fval' 'iter' 'kappa' 'kappa_opt' ...
    'label' 'm' 'max_Accu' 'maxiter' 'Mchol_1' 'Mchol_2' 'mean_Accu' 'mean_AUC' ...
    'max_report' 'mean_report' 'min_Accu' 'min_report' 'mj' 'mu' 'myoptions' ...
    'n' 'nu' 'output' 'perm' 'Prediction_X' 'report'  'seed' 'summary' 'theta'...
    'x' 'X' 'x0' 'xguess' 'Y'
    