%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB package for nonlinear Second-Order Cone programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,Gf,Hf]=f_svm_CoDo(x,theta)
% Function to minimize
% Support vector machine with Cobb-Douglas function
% f(w,b,ka_1,ka_2)=-theta*log(ka_1^2/(ka_1^2+1))-(1-theta)*log(ka_2^2/(ka_2^2+1))
% Input: 
%          x=(w,b,ka_1,ka_2) in R^{n+3}
%          theta := parameter of model
% Output:  f:= function 
%          Gf:= Gradient
%          HF:= Hessian

C1=theta;
C2=1-theta;
x=x(:);

n=length(x);
% Objective function
f=-C1*log(x(n-1)^2/(x(n-1)^2+1))-C2*log(x(n)^2/(x(n)^2+1));

% Compute the gradient
Gf=[zeros(n-3,1);0;-2*C1/((x(n-1)^2+1)*x(n-1));-2*C2/((x(n)^2+1)*x(n))];

% Hessian
z1=2*(3*x(n-1)^2+1)/((x(n-1)^2+1)*x(n-1))^2;
z2=2*(3*x(n)^2+1)/((x(n)^2+1)*x(n))^2;
Hf=sparse(n,n);
Hf(n-2,n-2)=0;
Hf(n-1,n-1)= C1*z1;
Hf(n,n)= C2*z2;

return
