%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDIPA : feasible directioninterior-point algorithm. See [1] for details
%
% The algorithm follows closely [1] page 1330.
% [1] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) A feasible 
%     direction algorithm for nonlinear second-order cone programs, 
%     Optimization Methods and Software, 34:6, 1322-1341, 
%     DOI: 10.1080/10556788.2018.1506452
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example1: Linear cone contraint Quadratic objective

% Construction of the matrices for the problem
N=16;
vertices=[ones(1,N); zeros(N-1,N)];
x0=10*eye(N);
x0=reshape(x0,N^2,1);
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
matQ=N*speye(N^2)-matJ;

my_options = fdipa_options('Display','iter','TolCon',1e-15,'Maxiter',300);
xmin = fdipaQuad(matQ,zeros(N^2,1),matA,b,x0,[],mj,my_options);

clear 'b'  'i' 'j' 'matJ' 'matA' 'matAi' 'matQ' 'mj' 'my_options' 'N' ...
    'vertices' 'x0' 'xmin' 'y0'