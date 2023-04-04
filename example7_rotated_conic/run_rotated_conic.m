%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDIPA : feasible direction interior-point algorithm
% An extension to NSOCP of the feasible direction interior-point algorithm 
% (FDIPA) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example: Quadratic objective with rotated conic constraint
seed = RandStream('mt19937ar','Seed',1);
m=10;
a = 10*rand(seed);
b = 10*rand(seed);
c = -1+ 2*rand(seed,m,1);
x0  = [a;b; rand(seed)*sqrt(2 *a*b)/norm(c)*c];

my_options = fdipa_options('Display','iter','TolCon',1e-15,'Maxiter',300);
xmin = fdipa(@(x)fun_rotated_conic(x,m),x0,[],@(x)g_rotated_conic(x,m),...
    [],my_options);

