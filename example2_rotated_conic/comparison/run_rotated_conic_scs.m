% Example 2: Quadratic objective with rotated conic constraint
% This example requires SCS solver
% https://www.cvxgrp.org/scs/index.html

%m=10;
%m=30;
%m=100;
m=1000;
%m=10000;
Q=eye(m+2);
Q(3,4)=-1;
Q(4,3)=-1;
data.P=2*sparse(Q);
A1=-[1 1 zeros(1,m); zeros(m,2) sqrt(2)*eye(m); 1 -1 zeros(1,m)];
data.A=sparse(A1);
data.b=zeros(m+2,1);
data.c=[1;1;2;1;zeros(m-2,1)];
cone.q=m+2;

% Optional solver settings
%settings = struct('eps_abs', 1e-9, 'eps_rel', 1e-9);

% Solve!
[x, y, s, info] = scs(data, cone);