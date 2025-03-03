% Example 2: Quadratic objective with rotated conic constraint
% This example requires CVX solver for matlab
% https://cvxr.com/cvx/doc/install.html

%m=10;
%m=30;
%m=100;
m=1000;
%m=10000;
A=[zeros(m,2) sqrt(2)*eye(m); 1 -1 zeros(1,m)];
Q=eye(m+2);
Q(3,4)=-1;
Q(4,3)=-1;

cvx_begin quiet
cvx_precision('low')
cvx_solver sedumi
%cvx_solver sdpt3
variable s(m+2) 
minimize (quad_form(s,Q)+s(1)+s(2)+2*s(3)+s(4))
subject to
norm(A*s)<=s(1)+s(2);

cvx_end



