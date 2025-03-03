% Example 1. Minimum distance between two ellipses as presented in 
% Examples 13.5 and 14.5 of [2]
% [2] A. Antoniou and W.-S. Lu. Practical Optimization: Algorithms and 
% Engineering Applications. Springer Publishing Company, Incorporated, 1st edition, 2007
% This example requires CVX solver for matlab
% https://cvxr.com/cvx/doc/install.html

A1=[0.5 0 0 0;0 1 0 0];
A2=[ 0 0 -0.7071 -0.7071; 0 0 -0.3536 0.3536];
b1=[-0.5;0];
b2=[4.2426;-0.7071];

cvx_begin quiet
cvx_precision('low')
%cvx_solver sdpt3
cvx_solver sedumi
variable x(4)
minimize ((x(1)-x(3))^2+(x(2)-x(4))^2)
subject to
norm(A1*x+b1)<=1;
norm(A2*x+b2)<=1;
cvx_end



