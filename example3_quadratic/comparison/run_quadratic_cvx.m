% Example 3: Linear cone contraint Quadratic objective
% This example requires CVX solver for matlab
% https://cvxr.com/cvx/doc/install.html

% Construction of the matrices for the problem
%N=16;
N=20;
%N=35;
%N=50;

vertices=[ones(1,N); zeros(N-1,N)];
b=-reshape(vertices,N^2,1);
matAi=eye(N);
matA=zeros(N^2);
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


cvx_begin quiet
cvx_precision('low')
cvx_solver sedumi
%cvx_solver sdpt3
variable x(N*N) 
minimize (0.5*quad_form(x,matQ))
subject to
for i=1:N
    A1=matA((i-1)*N+1,:);
    Ai=matA((i-1)*N+2:i*N,:);
    norm(Ai*x)<=A1*x+b((i-1)*N+1);
end
cvx_end



