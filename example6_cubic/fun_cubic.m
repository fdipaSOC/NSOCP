function [fun,grad_f] = fun_cubic(x)
% Example of cubic objective function
    x=x(:);
	fun = -x(1)*x(2)*x(3);
	grad_f = -[x(2)*x(3);x(1)*x(3);x(1)*x(2)];
