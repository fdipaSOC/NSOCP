%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,Gf,Hf]=f_svm_CoDo_L2(x,C1,C2)
% Function to minimize
% f(w,b,ka_1,ka_2)=0.5*|w|^2-C1*log(ka_1^2/(ka_1^2+1))-C2*log(ka_2^2/(ka_2^2+1))
% Input: 
%          x=(w,b,ka_1,ka_2) in R^{n+3}
%          C1, C2 := parameters of model
% Output:  f:= function 
%          Gf:= Gradient
%          HF:= Hessian

x=x(:);
n=length(x);
eps=1;
% Objective function
f=eps*0.5*norm(x(1:n-3))^2-C1*log(x(n-1)^2/(x(n-1)^2+1))-C2*log(x(n)^2/(x(n)^2+1));

% Compute the gradient
Gf=[eps*x(1:n-3);0;-2*C1/((x(n-1)^2+1)*x(n-1));-2*C2/((x(n)^2+1)*x(n))];

% Hessian
z1=2*(3*x(n-1)^2+1)/((x(n-1)^2+1)*x(n-1))^2;
z2=2*(3*x(n)^2+1)/((x(n)^2+1)*x(n))^2;
Hf=eps*speye(n,n);
Hf(n-2,n-2)=0;
Hf(n-1,n-1)= C1*z1;
Hf(n,n)= C2*z2;
return
