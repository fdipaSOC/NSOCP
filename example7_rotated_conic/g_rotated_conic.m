function [g_fun,grad_g]=g_rotated_conic(x,m)
% Rotated conic constraint using the change of variable described
% in ./documentation/doc_fdipa.pdf Example 7.7
    x=x(:);
    matA = zeros(m+2,m+2);
    matA(1,1) = 1;
    matA(1,2) = 1;
    matA(m+2,1) = 1;
    matA(m+2,2) = -1;
    matA(2:(m+1),3:m+2) = eye(m);
    g_fun = matA * x;
    grad_g = matA;