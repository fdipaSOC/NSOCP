%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDIPA : feasible directioninterior-point algorithm. See [1] for details
%
% The algorithm follows closely [1] page 1330.
% [1] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) A feasible 
%     direction algorithm for nonlinear second-order cone programs, 
%     Optimization Methods and Software, 34:6, 1322-1341, 
%     DOI: 10.1080/10556788.2018.1506452
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
my_options = fdipa_options('Display','iter');
[~,fval,~,output] =  fdipa(@fun_cubic,x0,@(x)gj_cubic(x,A,b),mj,y0,my_options);

%for paper [x,fval,exitflag,output]
fprintf('%d & %11f & %11.5e & %11f \\\\ \n', output.iterations,fval, output.firstorderopt, output.cputime)

clear 'A' 'b' 'x0' 'y0' 'mj' 'my_options'

function [fun,grad_f] = fun_cubic(x)
% Cubic objective function
    x=x(:);
    fun = -x(1)*x(2)*x(3);
    grad_f = -[x(2)*x(3);x(1)*x(3);x(1)*x(2)];
end
function [g_fun,grad_g] = gj_cubic(x,A,b)
% Linear constrain for the example of minimization of a cubic function
    x=x(:);
    g_fun = b + A*x;
    grad_g = A;
end
            
        