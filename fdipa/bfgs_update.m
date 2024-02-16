function hess = bfgs_update(xnew, x, grad_xnew, grad_x, hess)
% BFGS formula with the Han-Powell modification [1] Section 14.7
% INPUTS:
%   - xnew: point at which the updated hessian is needed.
%   - x: initial point at which the hessian is given.
%   - grad_xnew: gradient of the function at the point xnew.
%   - grad_x: gradient of the function at the point x.
%   - hess: hessian matrix at the point x.
% OUTPUTS:
%   - hess: approximate hessian matrix given by the BFGS method
% [1] D.G. Luenberger.Linear  and  Nonlinear  Programming:  Second  Edition.  Springer,2003.

	u = xnew-x; 
	v = grad_xnew-grad_x; 
	z = hess*u; 

    % theta parameter
	alpha = 0.2;
	if (u'*v >= alpha*u'*z)
		theta = 1;
	else
		theta = (1-alpha)*(u'*z)/(u'*z - u'*v);
	end
	r = theta*v + (1-theta)*z;
	gamma = (u'*r)/(u'*hess*u);
	hess = gamma*hess;

	% update step
	z = hess*u;
	hess = hess + (r*r'/(u'*r)) - (z*z'/(u'*z));

end



