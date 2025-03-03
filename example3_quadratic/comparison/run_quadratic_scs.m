% Example 3: Linear cone contraint Quadratic objective
% This example requires SCS solver
% https://www.cvxgrp.org/scs/index.html


% Construction of the matrices for the problem
%N=16;
N=20;
%N=35;
%N=50;

vertices=[ones(1,N); zeros(N-1,N)];
b=-reshape(vertices,N^2,1);
matAi=eye(N);
matA=zeros(N^2);
cvec=zeros(N^2,1);
for i=1:N
    matA((i-1)*N+1:i*N, (i-1)*N+1:i*N)=matAi(:,[2:i,1,(i+1):N]); 
end
mj=N*ones(1,N);
matJ=zeros(N^2);
for i=1:N
    for j=1:N
        matJ((i-1)*N+1:i*N,(j-1)*N+1:j*N)=speye(N);
    end
end
matQ=(N+1)*speye(N^2)-matJ;

data.P=sparse(matQ);
data.A=-sparse(matA);
data.b=b;
data.c=cvec;
cone.q=mj';

% Optional solver settings
%settings = struct('eps_abs', 1e-9, 'eps_rel', 1e-9);

% Solve!
[x, y, s, info] = scs(data, cone);