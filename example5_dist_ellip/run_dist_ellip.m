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

disp('Experiment 1: TolCon 1e-15, maximum iterations 300, bfgs Hessian update');
my_options = fdipa_options('Display','iter','TolCon',1e-15, ...
    'Maxiter',300,'Hessian','bfgs');
fdipa(@fun_dist_ellip,x0,y0,@g_dist_ellip,mj,my_options);

disp('Experiment 2: default options');
fdipa(@fun_dist_ellip,x0,y0,@g_dist_ellip,mj);

disp('Experiment 3: TolCon 1e-15, Hessian update off');
my_options = fdipa_options('Display','iter','TolCon',1e-15,'Hessian','off');
fdipa(@fun_dist_ellip,x0,y0,@g_dist_ellip,mj,my_options);

disp('Experiment 4: bfgs Hessian update, TolCon 1e-15');
my_options = fdipa_options('TolCon',1e-15,'Hessian','bfgs');
fdipa(@fun_dist_ellip,x0,y0,@g_dist_ellip,mj,my_options);

disp('Experiment 5: bfgs Hessian update');
my_options = fdipa_options('Hessian','bfgs');
fdipa(@fun_dist_ellip,x0,y0,@g_dist_ellip,mj,my_options);

clear 'x0' 'y0' 'mj' 'my_options'