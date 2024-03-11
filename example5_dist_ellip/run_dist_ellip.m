% In this example the problem of minimum distance between 
% two ellipses is considered
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

%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',2, output.iterations,fval, output.firstorderopt, output.cputime)

disp('Experiment 2:  Hessian update off');
my_options = fdipa_options('Display','final','HessianApproximation','off');
[~,fval,~,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',3, output.iterations,fval, output.firstorderopt, output.cputime)

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
                 