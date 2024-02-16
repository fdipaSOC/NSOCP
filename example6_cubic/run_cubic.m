% Example of cubic objective function with linear constraint
% min fun
% st  g_fun in K^1xK^1
% fun = -x(1)*x(2)*x(3)
% g_fun = b + A*x;

x0 = [10;10;10];
y0 = [];
mj = [1;1];

A = [1 2 2;-1 -2 -2];
b = [0;72];
    
my_options = fdipa_options('Display','iter','ConstraintTolerance',1e-15,'MaxIterations',50 ...
    ,'HessianApproximation','bfgs');
[~,fval,~,output] =  fdipa(@fun_cubic,x0,@(x)gj_cubic(x,A,b),mj,y0,my_options);

%for paper [x,fval,exitflag,output]
fprintf('%d & %11f & %11.5e & %11f \\\\ \n', output.iterations,fval, output.firstorderopt, output.cputime)

clear 'A' 'b' 'x0' 'y0' 'mj' 'my_options'