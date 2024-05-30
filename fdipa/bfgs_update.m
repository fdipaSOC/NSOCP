function hess = bfgs_update(x_new, x_old, y_new, y_old, fun, gj,hess_old)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	grad_x_new, grad_x_old
% BFGS formula with the Han-Powell modification [2] Section 14.7
% INPUTS:
%   - x_new: point at which the updated hessian is needed.
%   - x_old: initial point at which the hessian is given.
%   - y_new: lagrange multiplier at next iteration.
%   - y_old: lagrange multiplier at current iteration.
%   - fun:   handle to the objective function
%   - gj:    handle to the constraint function
%   - hess_old: hessian matrix at the point x.
% OUTPUTS:
%   - hess: approximate hessian matrix given by the BFGS method
% [2] D.G. Luenberger.Linear  and  Nonlinear  Programming:  Second  Edition.  Springer,2003.

    [~,grad_f_old] = fun(x_old); 
    grad_f_old = grad_f_old(:);
    [~,grad_g_old] = gj(x_old); 
    grad_x_old=grad_f_old -grad_g_old'*y_new;
    
	[~,grad_f_new] = fun(x_new);
    grad_f_new = grad_f_new(:);
    [~,grad_g_new] = gj(x_new);
    grad_x_new=grad_f_new -grad_g_new'*y_new;

	u = x_new-x_old; 
	v = grad_x_new-grad_x_old; 
	z = hess_old*u; 

    % theta parameter
	alpha = 0.2;
	if (u'*v >= alpha*u'*z)
		theta = 1;
	else
		theta = (1-alpha)*(u'*z)/(u'*z - u'*v);
	end
	r = theta*v + (1-theta)*z;
	gamma = (u'*r)/(u'*hess_old*u);
	hess = gamma*hess_old;

	% update step
	z = hess*u;
	hess = hess + (r*r'/(u'*r)) - (z*z'/(u'*z));

end



