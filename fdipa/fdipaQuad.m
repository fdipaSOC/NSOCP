function [x,fval,exitflag,Output,Lambda,grad,hessian] ...
    = fdipaQuad(Q,c,A,b,x0,mj,y0,myoptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This handler uses FDIPA to attempt to find a minimum of a 
% quadratic objective funcion with linear inequality constrains. 
% More specifically we solve for
% min 1/2 x' Q x + c x 
% s.t. A x  + b  in K^m1 x ... x K^mJ
% And the starting point for the optimization is x0.
% Note that fdipaQuad overrides the option for the hessian in such 
% a way the it enters as a fixed matrix to the algorithm by defining 
% the update hessian function as a constant.
%
% INPUTS:
%  - Q:  Symmetric positive definite matrix .
%  - c:  Vector with the coeffients of the linear term in the objective function.
%  - A:  Coefficient matrix of the inequality restrictions.
%  - b:  Right hand side of the inequality restrictions.
%  - x0: Initial condition.
%
% OPTIONAL INPUTS:
%  - mj: vector with the dimension of the cones associated with the 
%        restrictions, can be set to [] to use the default of the constraint 
%        that belong to a single cone.
%  - y0: initial feasible value for the lagrange multipliers, can be set to [].
%  - myoptions: Additional options passed to the solver. If empty default 
%        options are assumed. 
%
% OUTPUTS:
%  - x: argmin(fun)
%  - fval: Value of fun at the solution x.
%  - exitflag: Indicator of the status of the execution, has the value 0 
%              in successful.
%  - output:    - iterations:  Number of iterations taken
%               - cputime: Time used by the algorithm's execution.
%               - constrviolation:  A measure of how well the constraints are 
%               satisfied. A positive value indicates that all inequality 
%               constraints are satisfied strictly, and a negative value 
%               indicates that some of the constraints are violated by that 
%               amount.
%               - firstorderopt: Norm of the Lagrangian
%               - bestfeasible: Best solution found before the a stop 
%               condition is satisfied.
%               - stepsize: Length of last displacement in x.
%               - message: Exit message
%  - lambda: Structure with fields containing the Lagrange multipliers at 
%        the solution x
%  - grad: Gradient of Lagrangian at the solution x.
%  - hessian: Hessian of Lagrangian at the solution x.
%
% USAGE:
%   x = fdipaQuad(Q,c,A,b,x0,mj,y0,myoptions)
%   
%   [x,fval] = fdipaQuad(...)
%   [x,fval,exitflag] = fdipaQuad(...)
%   [x,fval,exitflag,output] = fdipaQuad(...)
%   [x,fval,exitflag,output,lambda] = fdipaQuad(...)
%   [x,fval,exitflag,output,lambda,grad] = fdipaQuad(...)
%   [x,fval,exitflag,output,lambda,grad,hessian] = fdipaQuad(...)
%
    %manage optional inputs
    if nargin >= 5
        if nargin < 6
            y0 = [];
        end
        if nargin < 7
            mj = length(x0);
        end
        if nargin < 8
            myoptions = fdipa_options();
        end
    else 
        error('fdipaQuad:Not enough input arguments')
    end
    % cast input vectors as column vectors
    x0 = x0(:);
    b = b(:); 
    c = c(:);
    if strcmp(myoptions.HessianApproximation,'default')
        myoptions.edit('HessianApproximation','user-supplied');
        b_update = @(x_new,x_old,y_new,y_old,fun,gj, hess_old) Q;
        myoptions.edit('HessianFcn',b_update);
        %check if the matrix is positive definite
        try chol(Q);
        catch 
            %if cholesky factorization fails indicates that matB is not positive definite.
            warning("The matrix Q is not positive definite, B^k = I will be used")
        end

    end
    
    [x,fval,exitflag,Output,Lambda,grad,hessian] ...
        =fdipa(@(x) quad_fun(x,Q,c),x0,@(x) linear_constraint_fun(x,A,b)...
        ,mj,y0,myoptions);    
end


function [f,grad_f]=quad_fun(x,Q,c)
    x = x(:);
    f=1/2 *x'*Q*x + c'*x;
    grad_f=Q*x + c;
end

function [g,grad_g]=linear_constraint_fun(x,A,b)
    x= x(:);
    g=A*x+b;
    grad_g=A;
end