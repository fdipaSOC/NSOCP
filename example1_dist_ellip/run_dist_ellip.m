%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example of fdipa
% Example 1. Minimum distance between two ellipses as presented in 
% Examples 13.5 and 14.5 of [2]
% [2] A. Antoniou and W.-S. Lu. Practical Optimization: Algorithms and 
% Engineering Applications. Springer Publishing Company, Incorporated, 1st edition, 2007
%
% $$\min f(x_1,x_2,x_3,x_4) $$
% Subject to
%  g1 in $K^3$, g2 in $K^3$
%
% $$f(x_1,x_2,x_3,x_4)=(x_1-x_3)^2+(x_2-x_4)^2$$
% $g1=[1;0.5*(x_1-1);x_2];$
% $g2=[1;-0.7071*x_3-0.7071*x_4+4.2426;-0.3536*x_3+0.3536*x_4-0.7071];$
%% Input Functions
% 
% First we need to define the objective function and its gradient 
% as a MATLAB(R) function, i.e.,
%
%    function [fun,grad_f]=fun_dist_ellip(x)
%    % Objective function for the example for distance between two ellipses
%    % Function to minimize
%    % f(x1,x2,x3,x4)=(x1-x3)^2+(x2-x4)^2
%        x=x(:);
%        fun =(x(1)-x(3))^2+(x(2)-x(4))^2;
%        grad_f=[2*(x(1)-x(3));2*(x(2)-x(4));-2*(x(1)-x(3));-2*(x(2)-x(4))];
%    end    
% 
% In equivalent way to the restrictions.
%
%    function [g_fun,grad_g]=g_dist_ellip(x)
%    % constraint function for the distance between two ellipses example
%        x= x(:);
%        g1=[1;0.5*(x(1)-1);x(2)];
%        g2=[1;-0.7071*x(3)-0.7071*x(4)+4.2426;...
%            -0.3536*x(3)+0.3536*x(4)-0.7071];
%        g_fun=[g1;g2];           
%        grad_g1=[0 0 0 0 ; 0.5 0 0 0; 0 1 0 0];
%        grad_g2=[0 0 0 0; 0 0 -0.7071 -0.7071;  0 0 -0.3536 0.3536];
%        grad_g=[grad_g1;grad_g2];
%    end                             
%
%% Options
% 
% To create and modify the options of the algorithm, we consider three 
% different experiments with different options. In the first one we use 
% all default options, so we can simply omit the parameter. For the 
% second one we want to only see a final summary of the algorithm, this 
% is done by setting |Display|  to |final| and for the Hessian approximation
% we want to disable the Hessian update (i.e. $B^k = I$), which 
% is done by setting the parameter |HessianApproximation| to |off|
%
%   myOptions = fdipa_options('Display','final','HessianApproximation','off');
%
% The last experiment is more involved, this time we want to only get the 
% a brief summary of the execution, so set |Display|  to |final|. Additionally
% we want to use a custom Hessian, ideally we would use the exact formula, 
% but since the objective function is not positive definite, that function 
% wont work as a matrix B^k, to solve this use add a multiple of the identity
% matrix to the exact Hessian, which guarantees that the resulting matrix is 
% positive definite. More precisely, if we want to use the matrix modnewton
%
%  mod_newton = [1 0 -1 0;0 1 0 -1;-1 0 1 0;0 -1 0 1];
%   epsilon = min(eig(mod_newton));
%   if epsilon <= 0.0
%     	mod_newton = mod_newton+(abs(epsilon)+0.1)*eye(4);
%   end
%
% as a custom Hessian, we have to set the option |HessianApproximation| to
% |user-supplied| and the option |HessianFcn| to be an anonymous function 
% with the appropriate inputs
%
%   @(x_new,x_old,y_new, y_old, fun,gj,hess_old) mod_newton
%
% Putting everything together we construct the options object using the following
%
%   my_options = fdipa_options('Display','final', ...
%       'HessianApproximation','user-supplied','HessianFcn',@(x_new,x_old,y_new, y_old, fun,gj,hess_old) mod_newton);
%
%% Other Inputs Arguments
%
% As we can see in this example we have 2 cone restrictions, whereby 
% $mj=\left[3;3\right]$.
% Finally you need to select a starting point |x0| and, if you want, 
% an approximation of the Lagrange multiplier |y0| . In this case we 
% are going to used the default value of |y0| .
% x0 = [1;0;1;5];

x0 = [1;0;1;5];
y0 = [];
mj = [3;3];
% if instead we want to use the subroutine searchStartingPoint to look for
% a feasible starting point
% x0=searchStartingPoint(4,@g_dist_ellip,mj);

disp('Experiment 1: default options');
[~,fval,~,output1] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0);
fprintf('BFGS & %d & %11f & %11f & %11.5e & %11f \\\\ \n', output1.iterations,fval, sqrt(fval), output1.firstorderopt, output1.walltime)

disp('Experiment 2:  Hessian update off');
my_options = fdipa_options('Display','final','HessianApproximation','off');
[~,fval,~,output2] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
fprintf('off & %d & %11f & %11f & %11.5e & %11f \\\\ \n', output2.iterations,fval, sqrt(fval), output2.firstorderopt, output2.walltime)

disp('Experiment 3:  Hessian update Modified Newton');
%for the modified newton approximation, because for this example the
%Hessian of the Lagrangian is constant, the Hessian approximation is also
%constant
mod_newton = [1 0 -1 0;0 1 0 -1;-1 0 1 0;0 -1 0 1];

try chol(mod_newton);
catch 
%if cholesky factorization fails indicates that matB is not positive definite.
    epsilon = min(eig(mod_newton));
  	mod_newton = mod_newton+(abs(epsilon)+0.1)*eye(4);
end

my_options = fdipa_options('Display','final', ...
    'HessianApproximation','user-supplied','HessianFcn',@(x_new,x_old,y_new, y_old, fun,gj,hess_old) mod_newton);
[~,fval,~,output3] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
fprintf('Modified Newton & %d & %11f & %11f & %11.5e & %11f \\\\ \n', output3.iterations,fval, sqrt(fval), output3.firstorderopt, output3.walltime)
%output for table 8 in the paper
fprintf('%d & %3.2f & %d & %3.2f & %d & %3.2f \n', output2.iterations, output2.walltime, output1.iterations, output1.walltime, output3.iterations, output3.walltime)


clear 'x0' 'y0' 'mj' 'my_options' 'x' 'fval' 'exitflag' 'output1' 'output2' 'output3' 'mod_newton' 'epsilon'

function [fun,grad_f]=fun_dist_ellip(x)
% Objective function for the example for distance between two ellipses
% Function to minimize
% f(x1,x2,x3,x4)=(x1-x3)^2+(x2-x4)^2
    x=x(:);
    fun =(x(1)-x(3))^2+(x(2)-x(4))^2;
    grad_f=[2*(x(1)-x(3));2*(x(2)-x(4));-2*(x(1)-x(3));-2*(x(2)-x(4))];
end    

function [g_fun,grad_g]=g_dist_ellip(x)
% constraint function for the distance between two ellipses example
    x= x(:);
    g1=[1;0.5*(x(1)-1);x(2)];
    g2=[1;-0.7071*x(3)-0.7071*x(4)+4.2426;...
        -0.3536*x(3)+0.3536*x(4)-0.7071];
    g_fun=[g1;g2];      
    
    grad_g1=[0 0 0 0 ; 0.5 0 0 0; 0 1 0 0];
    grad_g2=[0 0 0 0; 0 0 -0.7071 -0.7071;  0 0 -0.3536 0.3536];
    grad_g=[grad_g1;grad_g2];
end                             