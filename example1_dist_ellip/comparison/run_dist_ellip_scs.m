% Example 1. Minimum distance between two ellipses as presented in 
% Examples 13.5 and 14.5 of [2]
% [2] A. Antoniou and W.-S. Lu. Practical Optimization: Algorithms and 
% Engineering Applications. Springer Publishing Company, Incorporated, 1st edition, 2007
% This example requires SCS solver
% https://www.cvxgrp.org/scs/index.html

data.P = 2*sparse([1 0 -1 0;0 1 0 -1;-1 0 1 0;0 -1 0 1]);
A1=-[0 0 0 0;0.5 0 0 0;0 1 0 0];
A2=-[0 0 0 0; 0 0 -0.7071 -0.7071; 0 0 -0.3536 0.3536];
data.A=sparse([A1;A2]);
b1=[1;-0.5;0];
b2=[1;4.2426;-0.7071];
data.b=[b1;b2];
data.c=zeros(4,1);
cone.q=[3;3];

% Optional solver settings
%settings = struct('eps_abs', 1e-9, 'eps_rel', 1e-9);

% Solve!
[x, y, s, info] = scs(data, cone);




