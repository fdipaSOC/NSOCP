%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g,Gradg]=g_svm_ccl1(x,mu,Mchol_1,Mchol_2)
% Constraints for the SVM example written as connical restrictions
%
%
%  mu_1^t*w+b >= 1+ kappa_1* norm{S_1^t w}
%  -mu_2^t*w-b >= 1+  kappa_2* norm{S_2^t w}
%   kappa >=0
%  mu_1^t*w+b >= 0
%  -mu_2^t*w-b >= 0


x=x(:);
n=size(mu,2);

g1=[mu(1,:)*x(1:n)+x(n+1);x(n+2)*Mchol_1'*x(1:n)];
g2=[-mu(2,:)*x(1:n)-x(n+1);x(n+2)*Mchol_2'*x(1:n)];
g3=x(n+2);
g4=mu(1,:)*x(1:n)+x(n+1)-1;
g5=-mu(2,:)*x(1:n)-x(n+1)-1;
g=[g1;g2;g3;g4;g5];

% compute the gradient of constraints
Gradg1=[mu(1,:) 1 0;x(n+2)*Mchol_1' zeros(size(Mchol_1',1),1)  Mchol_1'*x(1:n)];
Gradg2=[-mu(2,:) -1 0;x(n+2)*Mchol_2' zeros(size(Mchol_2',1),1)  Mchol_2'*x(1:n)];
Gradg3=[zeros(1,n+1)  1];
Gradg4=[mu(1,:) 1 0];
Gradg5=[-mu(2,:) -1 0];
Gradg=[Gradg1;Gradg2;Gradg3;Gradg4;Gradg5];
return

     
