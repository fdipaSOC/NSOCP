function [fun,grad_f]=fun_rotated_conic(x,m)
% Quadratic objective function for the example with rotated conic constraint
    x= x(:);
    m = length(x)-2;
    fun= sum(x.^2)-2*x(3)* x(4) + x(1)+x(2)+2*x(3)+x(4);
    grad_f=2*x+[1;1;2-2*x(4);1-2*x(3);zeros(m-2,1)] ;