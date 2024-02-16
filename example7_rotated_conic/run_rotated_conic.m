%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDIPA : feasible direction interior-point algorithm
% An extension to NSOCP of the feasible direction interior-point algorithm 
% (FDIPA) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example: Quadratic objective with rotated conic constraint
seed = RandStream('mt19937ar','Seed',1);
m=1000;
a = 10*rand(seed);
b = 10*rand(seed);
c = -1+ 2*rand(seed,m,1);
x0  = [a;b; rand(seed)*sqrt(2 *a*b)/norm(c)*c];

my_options = fdipa_options('Display','iter', ...
    'ConstraintTolerance',1e-15,'MaxIterations',300);
[x,fval,exitflag,output] = fdipa(@(x)fun_rotated_conic(x,m),x0,@(x)g_rotated_conic(x,m),...
    [],[],my_options);

%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',m, output.iterations,fval, output.firstorderopt, output.cputime)