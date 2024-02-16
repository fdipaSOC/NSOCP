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
my_options = fdipa_options('Display','iter','ConstraintTolerance',1e-15, ...
    'MaxIterations',300,'HessianApproximation','bfgs');
[~,fval,~,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);

%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',1, output.iterations,fval, output.firstorderopt, output.cputime)

disp('Experiment 2: default options');
[~,fval,~,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0);

%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',2, output.iterations,fval, output.firstorderopt, output.cputime)

disp('Experiment 3: TolCon 1e-15, Hessian update off');
my_options = fdipa_options('Display','iter','ConstraintTolerance',1e-15,'HessianApproximation','off');
[~,fval,~,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',3, output.iterations,fval, output.firstorderopt, output.cputime)

disp('Experiment 4: bfgs Hessian update, TolCon 1e-15');
my_options = fdipa_options('ConstraintTolerance',1e-15,'HessianApproximation','bfgs');
[x,fval,exitflag,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',4, output.iterations,fval, output.firstorderopt, output.cputime)

disp('Experiment 5: bfgs Hessian update');
my_options = fdipa_options('HessianApproximation','bfgs');
[x,fval,exitflag,output] = fdipa(@fun_dist_ellip,x0,@g_dist_ellip,mj,y0,my_options);
%for paper [x,fval,exitflag,output]
fprintf('%d & %d & %11f & %11.5e & %11f \\\\ \n',5, output.iterations,fval, output.firstorderopt, output.cputime)

clear 'x0' 'y0' 'mj' 'my_options' 'x' 'fval' 'exitflag' 'output'