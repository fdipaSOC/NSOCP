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

     