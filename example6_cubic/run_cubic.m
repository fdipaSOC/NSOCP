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
    
my_options = fdipa_options('Display','iter','TolCon',1e-15,'Maxiter',50 ...
    ,'Hessian','bfgs');
xmin=fdipa(@fun_cubic,x0,y0,@(x)gj_cubic(x,A,b),mj,my_options);

clear 'A' 'b' 'x0' 'y0' 'mj' 'xmin' 'my_options'