%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1. Minimum distance between two ellipses as presented in 
% Examples 13.5 and 14.5 of [2]
% [2] A. Antoniou and W.-S. Lu. Practical Optimization: Algorithms and 
% Engineering Applications. Springer Publishing Company, Incorporated, 1st edition, 2007
%
% min f(x1,x2,x3,x4) 
% st g1 in K^3, g2 in K^3
%
% f(x1,x2,x3,x4)=(x1-x3)^2+(x2-x4)^2
% g1=[1;0.5*(x(1)-1);x(2)];
% g2=[1;-0.7071*x(3)-0.7071*x(4)+4.2426;-0.3536*x(3)+0.3536*x(4)-0.7071];

x0 = [1;0;1;5];
y0 = [];
mj = [3;3];

disp('Experiment 1: default options');
[~,fval,~,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0);
fprintf(' BFGS & %d & %11f & %11f & %11.5e & %11f \\\\ \n', output.iterations,fval, sqrt(fval), output.firstorderopt, output.walltime)

disp('Experiment 2:  Hessian update off');
my_options = fdipa_options('Display','final','HessianApproximation','off');
[~,fval,~,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
fprintf(' off & %d & %11f & %11f & %11.5e & %11f \\\\ \n', output.iterations,fval, sqrt(fval), output.firstorderopt, output.walltime)

disp('Experiment 3:  Hessian update Modified Newton');
%for the modified newton approximation, because for this example the
%Hessian of the Lagrangian is constant, the Hessian approximation is also
%constant
mod_newton = [1 0 -1 0;0 1 0 -1;-1 0 1 0;0 -1 0 1];
epsilon = min(eig(mod_newton));
if epsilon <= 0.0
  	mod_newton = mod_newton+(abs(epsilon)+0.1)*eye(4);
end

my_options = fdipa_options('Display','final', ...
    'HessianApproximation','user-supplied','HessianFcn',@(x_new,x_old,y_new, y_old, fun,gj,hess_old) mod_newton);
[~,fval,~,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
fprintf(' Modified Newton & %d & %11f & %11f & %11.5e & %11f \\\\ \n', output.iterations,fval, sqrt(fval), output.firstorderopt, output.walltime)


clear 'x0' 'y0' 'mj' 'my_options' 'x' 'fval' 'exitflag' 'output'

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