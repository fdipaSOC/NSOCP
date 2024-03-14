% Example: Quadratic objective with rotated conic constraint
%
% FDIPA : feasible directioninterior-point algorithm. See [1] for details
%
% The algorithm follows closely [1] page 1330.
% [1] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) A feasible 
%     direction algorithm for nonlinear second-order cone programs, 
%     Optimization Methods and Software, 34:6, 1322-1341, 
%     DOI: 10.1080/10556788.2018.1506452


seed = RandStream('mt19937ar','Seed',1);
m=1000      0;
a = 10*rand(seed);
b = 10*rand(seed);
c = -1+ 2*rand(seed,m,1);
x0  = [a;b; rand(seed)*sqrt(2 *a*b)/norm(c)*c];

my_options = fdipa_options('Display','final');
[~,fval,~,output] = fdipa(@(x)fun_rotated_conic(x,m),x0,@(x)g_rotated_conic(x,m),...
    [],[],my_options);

%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',m, output.iterations,fval, output.firstorderopt, output.cputime)

clear 'seed' 'a' 'b' 'c' 'x0' 'myoptions' 'fval' 'output'

function [fun,grad_f]=fun_rotated_conic(x,m)
% Quadratic objective function for the example with rotated conic constraint
    x= x(:);
    m = length(x)-2;
    fun= sum(x.^2)-2*x(3)* x(4) + x(1)+x(2)+2*x(3)+x(4);
    grad_f=2*x+[1;1;2-2*x(4);1-2*x(3);zeros(m-2,1)] ;
end

function [g_fun,grad_g]=g_rotated_conic(x,m)
% contraint function for the rotated conic constraint
    x=x(:);
    matA = zeros(m+2,m+2);
    matA(1,1) = 1;
    matA(1,2) = 1;
    matA(m+2,1) = 1;
    matA(m+2,2) = -1;
    matA(2:(m+1),3:m+2) = eye(m);
    g_fun = matA * x;
    grad_g = matA;
end