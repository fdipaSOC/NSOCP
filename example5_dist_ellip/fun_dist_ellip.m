function [fun,grad_f]=fun_dist_ellip(x)
% Objective function for the example for distance between two ellipses
% Function to minimize
% f(x1,x2,x3,x4)=(x1-x3)^2+(x2-x4)^2
    x=x(:);
    fun =(x(1)-x(3))^2+(x(2)-x(4))^2;
    grad_f=[2*(x(1)-x(3));2*(x(2)-x(4));-2*(x(1)-x(3));-2*(x(2)-x(4))];

    