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


%%experiment 1 : logarithmic objective with weights 1/2, 1/2
disp('Experiment 1: logarithmic objective with weights 1/2, 1/2, BFGS');
% use BFGS hessian update for the optimzation
myoptions = fdipa_options('MaxIterations',10000,'Display','final',...
    'NumericalConditioning',1e15, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10);


fun=@(x)f_svm_CoDo(x,theta);
const=@(x)g_svm(x,mu,Mchol_1,Mchol_2);

for iter=1:maxiter
    nu = 0.1*rand(seed,2,1);
    kappa=sqrt(nu./(1-nu));
    
    %look for feasible starting point using an auxiliary problem
    %in order to find different feasible points, we imporse 
    % lower bounds in the values of kappa which are selecte randomly
    x0_aux= [zeros(n+3,1);10];
    myoptions_aux = fdipa_options('Display','off','LowerOptimalityBound', 0);
    x_aux=fdipa(@(z) f_svm_aux(z),x0_aux,@(z) g_svm_aux(z,mj,mu,kappa,Mchol_1,Mchol_2),mj,[],myoptions_aux);
    x0 = x_aux(1:(end-1));

    % We search for a feasible starting point using CVX
    % cvx_begin quiet   
    % variables w(n) b
    % subject to
    % kappa*norm(Mchol_1'*w,2)<=w'*mu(1,:)'+b-1;
    % kappa*norm(Mchol_2'*w,2)<=-w'*mu(2,:)'-b-1;
    % cvx_end
    % x0=[w;b;kappa];
    
    [x,fval,~,output]=fdipa(fun,x0,const,mj,[],myoptions);
    kappa_opt=x(n+2:end);
    eta_opt=kappa_opt.^2./(1+kappa_opt.^2);    
    report(iter,:) = [output.iterations, output.walltime, fval,output.firstorderopt , output.compslack,output.constrviolation,eta_opt(1),eta_opt(2)];
    fprintf('%d & %d & %11f & %11.5e & %11f & %11.5f & %11.5f \\\\ \n',m, output.iterations,fval, output.firstorderopt, output.walltime,eta_opt(1),eta_opt(2))
   
end
max_report = max(report);
min_report = min(report);
mean_report = mean(report);
array2table([max_report; min_report;mean_report] ,...
    'VariableNames',{'iterations','time','fval','norm_lag','comp_slack','feasibility','eta1', 'eta2'},...
    'RowNames',{'max','min','mean'})

fprintf("{\\bf " + label + "} & %4.1f & %d & %d & %4.2f & %4.2f & %4.2f & %4.4f & %4.4f & %11.5e & %11.5e \\\\\n", ...
mean_report(1), min_report(1), max_report(1),mean_report(2), min_report(2), max_report(2),report(1,7),report(1,8),report(1,4),max_report(1,5));

function [f,Gf]=f_svm_aux(x)
    f = x(end);
    Gf = [zeros(length(x)-1,1);1];
end

function [g,Gradg]=g_svm_aux(x,mj,mu,kappa,Mchol_1,Mchol_2)
    [g,Gradg] = g_svm(x(1:(end-1)),mu,Mchol_1,Mchol_2);
    %impose a lower bound in the values of the kappa vector
    g(end-1:end)= g(end-1:end)-kappa;
    Gradg= [Gradg,zeros(sum(mj),1)];

    n_cones = length(mj); % number of cones
    block_begin = ones(n_cones,1); % index of the first coordinate of i-th cone
    block_end = mj; % index of the last coordinate of i-th cone 
    if n_cones>1
        for i=2:n_cones
            block_end(i)=block_end(i-1)+mj(i);
            block_begin(i)=block_end(i-1)+1;
        end
    end 

    for k=1:n_cones
        g(block_begin(k)) = g(block_begin(k)) + x(end);
        Gradg(block_begin(k),:) = Gradg(block_begin(k),:)+ [zeros(1,length(x)-1), 1];
    end
end


