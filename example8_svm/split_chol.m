%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB package for nonlinear Second-Order Cone programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu,Mchol_1,Mchol_2]=split_chol(X,Y)
% Separates the data using the given classifications,
% Measure something and perform a Cholesky factorization
Min_label=min(Y);
if Min_label<0
    find1=find(Y==1);
    find2=find(Y==-1);
else
    find1=find(Y==1);
    find2=find(Y==2);
end
n=size(X,2);
x=X(find1,:);
xx=X(find2,:);
mu(1,:)=mean(x);
mu(2,:)=mean(xx);
Sigma(:,:,1)=cov(x)+1e-6*eye(n);
Sigma(:,:,2)=cov(xx)+1e-6*eye(n);
Mchol_1=chol(Sigma(:,:,1),'lower');
Mchol_2=chol(Sigma(:,:,2),'lower');