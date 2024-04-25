%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDIPA : feasible direction interior-point algorithm. 
% The details of this particular implementation are contained in [1]. 
% The algorithm was first presented in [2, page 1330]. This is an extension 
% to NSOCP of the feasible direction interior-point algorithm (FDIPA) 
% proposed by Herskovits in [2] for smooth constrained nonlinear programming (NP).
%
% This package has been download from https://github.com/fdipaSOC/NSOCP
% See README.md for details.
% 
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB package for nonlinear Second-Order Cone programs
% [2] Alfredo Canelas, Miguel Carrasco, Julio Lopez (2019) A feasible 
%     direction algorithm for nonlinear second-order cone programs, 
%     Optimization Methods and Software, 34:6, 1322-1341, 
%     DOI: 10.1080/10556788.2018.1506452
% [3] J. Herskovits, Feasible direction interior-point technique for 
%     nonlinear optimization, J. Optim. TheoryAppl. 99 (1998), pp. 121–146.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,fval,exitflag,output,lambda,grad,hessian] = fdipa(varargin)
% FDIPA : feasible direction interior-point algorithm
% fdipa attempts to find a minimum of a scalar function with 
% second order cone constraints, starting at an initial estimate.
%
% INPUTS:
%  - fun: Objective function, the gradient must be provided as a second output.
%  - x0:  Starting point for the optimization.
%  - gj:  The function that computes the cone constraint gj. gj accepts a 
%         vector x and returns a vector gj with the evaluation of the constraint 
%         function as first output and the gradient of the constraints as 
%         second output. See documentation for details.
%
% OPTIONAL INPUTS:
%  - mj:  If the problem has more than one cone constraint, then 
%         gj=[g1;...;gJ] and the dimension of each cone should be specified 
%         through the vector mj=[m1,...,mJ].  
%  - y0:  An initial estimate of the lagrange multiplier assiciated with the 
%         cone constraints. This argument is optional, if you do not want to 
%         provide it set y0 = [], in this case the default value of y0 is 
%         (1, 0, 0, ... , 0) for each cone.
%  - options: Object created by fdipa_options() containing the configurations 
%         for the execution of the algorithm.
%
% ALTERNATIVE INPUTS:
%  - problem: Structure containing the information for the execution of the 
%          algorithm.
%        problem.objective: Objective function
%        problem.x0: Initial point for x   
%        problem.gj: The funcion that computes the cone constraint gj. 
%        
%   OPTIONAL:
%        problem.mj: Vector with the dimensions of the cones constraints 
%        problem.y0: Initial point for the Lagrange multiplier
%        problem.options: Options created with fdipa_options()
%
% OUTPUTS:
%  - x: argmin(fun)
%  - fval: Value of fun at the solution x.
%  - exitflag: Indicator of the status of the execution, has the value 0 
%              in successful.
%  - output: A structure with information about the optimization. The 
%            fields of the structure are the following:
%        output.iterations: Number of iterations taken
%        output.constrviolation: A nonnegative number that measures 
%            the feasibility of the solution. It is identically zero for a 
%            feasible solution.
%        output.stepsize: Length of last displacement in x
%        output.firstorderopt: Measure of first-order optimality given by the 
%            norm of the gradient of the Lagrangian
%        output.message: Exit message
%        output.compslack: A measure of the complementary slackness in the 
%            KKT conditions
%        output.bestfeasible: structure with the best feasible solution
%            found
%  - lambda: Structure with fields containing the Lagrange multipliers at the 
%        solution x:
%        lambda.gj:  lagrange multiplier corresponding to the inequality contraints
%  - grad: Gradient of fun at the solution x.
%  - hessian: Hessian approximation at the solution x.
%
% USAGE:
%   x = fdipa(fun,x0,gj)
%   x = fdipa(fun,x0,gj,mj)
%   x = fdipa(fun,x0,gj,mj,y0)
%   x = fdipa(fun,x0,gj,mj,y0,options)
%   x = fdipa(problem)
%   
%   [x,fval] = fdipa(...)
%   [x,fval,exitflag] = fdipa(...)
%   [x,fval,exitflag,output] = fdipa(...)
%   [x,fval,exitflag,output,lambda] = fdipa(...)
%   [x,fval,exitflag,output,lambda,grad] = fdipa(...)
%   [x,fval,exitflag,output,lambda,grad,hessian] = fdipa(...)

    tic
    exitflag = 0;
    % verification for the number of arguments is correct
    if ~any([1,3,4,5,6]==nargin)
        error('fdipa: Wrong number of inputs.');
    end
    % Sort the input data in the corresponding variables
    if nargin==1
        if ~isstruct(varargin{1})
            error('fdipa: The input has not the correct format');
        end
        problem = varargin{1};
        if ~all(isfield(problem,{'objective','x0','gj'}))
            error('fdipa: The input structure is missing required fields');
        end
        fun = problem.objective;
        x0 = problem.x0;
        gj = problem.gj;
        if isfield(problem,'mj')
            mj = problem.mj;
        else
            mj = [];
        end
        if isfield(problem,'y0')
            y0 = problem.y0;
        else
            y0 = [];
        end
        if isfield(problem,'options')
            options = problem.options;
        else
            options = options_class();
        end
    end

    if nargin >= 3
        fun = varargin{1};
        x0 = varargin{2};
        gj = varargin{3};
        if nargin >= 4
            mj = varargin{4};
            if nargin >=5
                y0=varargin{5};
                if nargin >=6
                    options = varargin{6};
                else
                    options = options_class();
                end
            else
                y0=[];
                options = options_class();
            end
        else
            mj=[];
            y0=[];
            options = options_class();
        end
    end
  
    % preliminary checks on the input data
    % check if the gradient of the objective function is provided
    if nargout(fun)==1
        error('fdipa: You must supply the gradient of the objetive function.');
    end
    
    % check if the gradient of the constraints gj provided
    if nargout(gj)==1
        error('fdipa: You must supply the gradient of the cone constraints.');
    end
    
    % set the initial condition to a column vector
    x0=x0(:);
    
    % check compatibility of dimensions between the input variables and functions
    [~,grad_fx0] = fun(x0);   
    [gx0,grad_gx0] = gj(x0); 
    gx0 = gx0(:);
    dimx = length(x0);   
    if length(grad_fx0)~=dimx
        error('fdipa: grad_fx has incompatible dimensions.');
    end     
    
    % if the dimension of the cones "mj" are not provided, it is 
    % assumed that there is only one cone
    if min(size(mj))==0
        mj = length(gx0);
    end
    %the dimensions of the cones must agree with the contraints
    if length(gx0) ~= sum(mj)
        error('fdipa: The dimensions of the cones mj are incompatible with gj(x).');
    end
    
    n_cones = length(mj); % number of cones
    block_begin = ones(n_cones,1); % index of the first coordinate of i-th cone
    block_end = mj; % index of the last coordinate of i-th cone 
    if n_cones>1
        for i=2:n_cones
            block_end(i)=block_end(i-1)+mj(i);
            block_begin(i)=block_end(i-1)+1;
        end
    end 
    
    % check if y0 is provided. By default we take the direction (1,0,...,0) on 
    % each cone, which is always feasible
    if min(size(y0))==0
        y0 = zeros(length(gx0),1);
        for i=1:n_cones
            y0(block_begin(i))=1;
        end
    end
    dimg = length(y0); 
      
    %check if the contraint function has the correct dimensions
    if length(gx0)~=dimg
        error('fdipa: gj(x) has incompatible dimensions.');
    end 
    
    [nrows_gradg, ncol_gradg]=size(grad_gx0);
    if ncol_gradg~=dimx | nrows_gradg~=dimg
        error('fdipa: grad_gj has incompatible dimensions.');
    end 

    % End of validations    
    
    % Options load: Hessian
    if strcmp(options.HessianApproximation,'default')
        options.HessianApproximation = 'bfgs';
    end

    switch options.HessianApproximation
        case 'off'
            b_update = @(x_new,x_old,y_new,y_old,fun,gj,hess_old) eye(dimx);
        case 'bfgs'
            b_update = @bfgs_update;
            if isnan(options.HessianResetIterations)
                options.HessianResetIterations = dimx;
            end
        case 'user-supplied'
            if ~isa(options.HessianFcn,'function_handle')
                error(append('fdipa: options.Hessian = supplied, you ', ...
                    'must supply an update of the Hessian or change ', ...
                    'options.HessianApproximation'));
            else
                b_update = options.HessianFcn;
            end
        otherwise 
            error('fdipa:Choosen option for the Hessian is not valid.');
    end

    if isnan(options.HessianResetIterations)
        options.HessianResetIterations = options.MaxIterations+1 ;
    end

    % Start of the algorithm    
    % Step 0
    xk = x0; 
    yk = y0;
    matB=eye(dimx); 
    lamb_min =spectral_decomposition(gx0,mj);
    % Check feasible point, min(lamb_min) < 0 means that the point is unfeasible
    %if min(lamb_min)<-options.ConstraintTolerance(1)
    if min(lamb_min)<0
        x = x0;
        [fval,grad] = fun(x); 
        exitflag = -1;
        % output
        output.iterations = 0;
        output.cputime = toc;
        output.constrviolation = max(abs(min(lamb_min, 0)));
        output.stepsize = [];
        output.firstorderopt = [];
        output.message = 'fdipa:starting point is unfeasible';
        output.bestfeasible=[];
        % lambda
        lambda.gj = y0;
        % Hessian
        hessian = matB;
        disp(output.message);
        return
    end

    for k=0:options.MaxIterations
        % Step 1
        % evaluate Lagrangian at (x,yk)
        [gxk,grad_gxk] = gj(xk); 
        gxk = gxk(:);       
        [fxk,grad_fxk] = fun(xk);
        grad_fxk = grad_fxk(:);
        norm_lag = norm(grad_fxk -grad_gxk'*yk);  

        %Report iteration information
        % Start screen information table in case 'iter'             
        if strcmp(options.Display,'iter')
            % compute complementary slack
            slack_cones =zeros(n_cones,1);
            for i=1:n_cones
                slack_cones(i) = abs(dot(gxk(block_begin(i):block_end(i)) ...
                                             ,yk(block_begin(i):block_end(i))));
            end
            comp_slack = max(slack_cones);
            % verify feasibility
            if k==0
                lamb_min_t = spectral_decomposition(gxk,mj);
            end
            feasibility=max(abs(min(lamb_min_t, 0)));

            if k==0
                fprintf(append(repelem(' ',32),'Feasibility  First-order    Dual step  Comp. Slack.\n'));
                fprintf(append('Iter         f(x)    Step-size                optimality  ',...
                    '  |yk+1-yk|   max|gj.yj|\n'));
                fprintf(['%-4d  %11.6e  ',repelem(' ',13),'%11.3e  %11.3e  ',repelem(' ',13),'%11.3e\n'],...
                    k,fun(xk), feasibility, norm_lag, comp_slack);
            else
                fprintf('%-4d  %11.6e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e\n', ...  
                    k,fun(xnew),stepsize,feasibility,...
                    norm_lag,dual_stepsize, comp_slack);
            end
        end
        if k == options.MaxIterations
            break
        end
   
        % stop if first order optimallity is satified
        if norm_lag <options.OptimalityTolerance
            break 
        end
        if k>0 && stepsize < options.StepTolerance
            break
        end

        % Step 1: looking for a feasible descend direction. 
        % Solving systems (i) and (ii)
        % depending on the conditioning of the systems a different approach 
        % for solving the linear system might be more accurate
        % if the system is too ill-conditioned, the Hessian approximation 
        % will be reset to the identity before continuing
        arrw_y=arrow(yk,mj); % Arrow matrix of yk
        arrw_g=arrow(gxk,mj); % Arrow matrix of gj(xk)      

        %is_linear_system_solved = false;
        % Solving system for (da, ya)
        [daya,rcond1]=linsolve([matB -grad_gxk'; arrw_y*grad_gxk arrw_g],[-grad_fxk; zeros(sum(mj),1)]);
        da = daya(1:dimx);
        ya = daya((dimx+1):(dimx + sum(mj)));
        % Solving system for (db, yb)
        [dbyb,rcond2] = linsolve([matB,-grad_gxk';arrw_y*grad_gxk, arrw_g],[zeros(dimx,1);yk]);
        db = dbyb(1:dimx); 
        yb = dbyb((dimx+1):end);    
        
        %check if solution of the linear systems is accurate
        if rcond1 < 1/options.NumericalConditioning || rcond2 < 1/options.NumericalConditioning
            is_linear_system_solved =  norm([matB, -grad_gxk'; arrw_y*grad_gxk, arrw_g]*[da;ya] - [-grad_fxk; zeros(sum(mj),1)]) < options.LinearSystemTolerance ...
                && norm([matB,-grad_gxk';arrw_y*grad_gxk, arrw_g]*[db;yb] - [zeros(dimx,1);yk]) < options.LinearSystemTolerance;
            %norm([matB, -grad_gxk'; arrw_y*grad_gxk, arrw_g]*[da;ya] - [-grad_fxk; zeros(sum(mj),1)])
            %norm([matB,-grad_gxk';arrw_y*grad_gxk, arrw_g]*[db;yb] - [zeros(dimx,1);yk])
        else
            is_linear_system_solved = true;
        end        
        % if the first approach fail, we use smaller matrices by solving 
        % one of the equations manually
        if ~is_linear_system_solved
            % Solving system for (da, ya)
            arrw_inv_g = arrow_inv(gxk,mj);
            matM1=arrw_inv_g*arrw_y*grad_gxk;
            [da,rcond1]=linsolve((matB+grad_gxk'*matM1),-grad_fxk);
            ya=-matM1*da;
            % Solving system for (db, yb)
            matM2 = arrw_y * grad_gxk;
            [db,rcond2] = linsolve((matB - grad_gxk' * arrw_inv_g * matM2),(grad_gxk' * arrw_inv_g * yk));
            yb = arrw_inv_g *(yk - matM2 * db );
            %check if solution of the linear systems is accurate
            if rcond1 < 1/options.NumericalConditioning || rcond2 < 1/options.NumericalConditioning
                is_linear_system_solved =  norm([matB, -grad_gxk'; arrw_y*grad_gxk, arrw_g]*[da;ya] - [-grad_fxk; zeros(sum(mj),1)]) < options.LinearSystemTolerance ...
                    && norm([matB,-grad_gxk';arrw_y*grad_gxk, arrw_g]*[db;yb] - [zeros(dimx,1);yk]) < options.LinearSystemTolerance;
                %norm([matB, -grad_gxk'; arrw_y*grad_gxk, arrw_g]*[da;ya] - [-grad_fxk; zeros(sum(mj),1)])
                %norm([matB,-grad_gxk';arrw_y*grad_gxk, arrw_g]*[db;yb] - [zeros(dimx,1);yk])
            else
                is_linear_system_solved = true;
            end 
        end
        % if the second approach fail, it means that the conditioning of 
        % the system is too bad for the current choice of the Hessian 
        % approximation matB, so we reset it to matB = Id and try again.
        if ~is_linear_system_solved
            matB = eye(dimx);
            % Solving system for (da, ya)
            arrw_inv_g = arrow_inv(gxk,mj);
            matM1=arrw_inv_g*arrw_y*grad_gxk;
            [da,rcond1]=linsolve((matB+grad_gxk'*matM1),-grad_fxk);
            ya=-matM1*da;
            % Solving system for (db, yb)
            matM2 = arrw_y * grad_gxk;
            [db,rcond2] = linsolve((matB - grad_gxk' * arrw_inv_g * matM2),(grad_gxk' * arrw_inv_g * yk));
            yb = arrw_inv_g *(yk - matM2 * db );
            %check if solution of the linear systems is accurate
            if rcond1 < 1/options.NumericalConditioning || rcond2 < 1/options.NumericalConditioning
                is_linear_system_solved =  norm([matB, -grad_gxk'; arrw_y*grad_gxk, arrw_g]*[da;ya] - [-grad_fxk; zeros(sum(mj),1)]) < options.LinearSystemTolerance ...
                    && norm([matB,-grad_gxk';arrw_y*grad_gxk, arrw_g]*[db;yb] - [zeros(dimx,1);yk]) < options.LinearSystemTolerance;
                norm([matB, -grad_gxk'; arrw_y*grad_gxk, arrw_g]*[da;ya] - [-grad_fxk; zeros(sum(mj),1)])
                norm([matB,-grad_gxk';arrw_y*grad_gxk, arrw_g]*[db;yb] - [zeros(dimx,1);yk])
            else
                is_linear_system_solved = true;
            end            
        end
        % If we cannot solve the linear system accurately the algorithm cannot continue
        if ~is_linear_system_solved
            break
        end
        % If norm(da) = 0 we are at a KKT point
        norm_da = norm(da); 
        if norm_da < options.StepTolerance
            break
        end        

        % (iii) Compute rho        
        denom_rho = db'*grad_fxk;
        if denom_rho>0
            rho = min(options.ParPhi*norm_da^2,...
                 (options.ParXi-1)*(da'*grad_fxk)/denom_rho); 
        else
            rho = options.ParPhi*norm_da^2; 
        end

        % (iv) Compute feasible descent direction
        desc_dir = da + rho*db;  
        yk_hat = ya + rho*yb;
        
        
        % Step 2: Armijo linear search
        % to accept a step size we check if the new point has enough descend, 
        % the point we reach is feasible, and an additional 
        % technical spectral condition specified in the algorithm 
        % which is related to the constraints      
        
        t=1;
        [lambda1_gxk,lambda2_gxk,u1_gxk,u2_gxk] ...
            = spectral_decomposition(gxk,mj);

        first_armijo_iteration = true;
        % armijo_iter = 1;
        while first_armijo_iteration || (~has_enough_descend)||(~is_feasible) ...
                ||(~satisfy_spectral_condition)
            % armijo_iter = armijo_iter + 1;
            if first_armijo_iteration
                first_armijo_iteration = false;
            else
                t=t*options.ParNu; 
            end
            fxt=fun(xk + t*desc_dir);
            gxt=gj(xk + t*desc_dir); 
            lamb_min_t=spectral_decomposition(gxt,mj);
            has_enough_descend = fxt<=fxk+t*options.ParEta*desc_dir'*grad_fxk;
            %is_feasible = min(lamb_min_t)>-options.ConstraintTolerance(1);
            is_feasible = min(lamb_min_t)>0;
            % verify spectral conditions
            satisfy_spectral_condition = true;
            [lambda1_gxt,lambda2_gxt] ...
                = spectral_decomposition(gxt,mj);
            %This is done for each cone
            for j = 1:n_cones
                range = block_begin(j):block_end(j);
                if (u1_gxk(range)'*yk_hat(range) < 0) ...
                        && (lambda1_gxk(j) < options.ParLambdam(1) )
                    if lambda1_gxt(j) < lambda1_gxk(j)
                        satisfy_spectral_condition = false;
                        break
                    end
                end
                
                if (u2_gxk(range)'*yk_hat(range) < 0) ...
                        && (lambda2_gxk(j) < options.ParLambdam(1) )
                    if lambda2_gxt(j) < lambda2_gxk(j) 
                        satisfy_spectral_condition = false;
                        break
                    end
                end
            end
        end
        
        % Step 3: Updates
        % New primal point
        xnew=xk + t*desc_dir;       
        % update y 
        ynew = zeros(dimg,1);
        [~,~,u1_gxt,u2_gxt] = spectral_decomposition(gxt,mj);
        for i=1:n_cones
            if mj(i)==1
                ynew(block_begin(i))= max(min(options.ParCS,ya(block_begin(i))),options.ParCI);
            else
                % Projection of yak in the linear subspace generated by u1, u2
                % tilde_yak = yu1 * u1 + yu2 * u2
                % next project tilde_yak in convex region of admissible values 
                % for the Lagrange multiplier. We multiply by 2 because 
                % u1,u2 are not unit vectors. By construction the vector ynew 
                % always lies in the cone of admissible values for the multiplier.
                yu1 = 2*dot(ya(block_begin(i):block_end(i)) ...
                        ,u1_gxt(block_begin(i):block_end(i)));
                yu2 = 2*dot(ya(block_begin(i):block_end(i)) ...
                        ,u2_gxt(block_begin(i):block_end(i)));
                % update y
                ynew(block_begin(i):block_end(i)) ...
                    = max(options.ParCI,min(options.ParCS,yu1)) ...
                        *u1_gxt(block_begin(i):block_end(i)) ...
                    +max(options.ParCI,min(options.ParCS,yu2)) ...
                        *u2_gxt(block_begin(i):block_end(i));
            end    
        end
        
        % Update Bk
        matB = b_update(xnew,xk,ynew,yk,fun,gj, matB);

        % To verify Assumtion 3.5 in [1] we compute the largest singular 
        % value norm(matB,2) and the condition number cond(matB) of the 
        % matrix matB, compare them with the parameters options.ParSigma1 
        % and options.ParSigma2. This is a measure of the conditioning 
        % of the matrix. If it fails the test the matrix B to the identity. 

        % Condition on smallest eigenvalue not being too small
        if norm(matB,2)/cond(matB) < options.ParSigma1 
            matB = eye(dimx);
        end
        % Condition on largest singular value not being too large
        if norm(matB,2) > options.ParSigma2 
            matB = eye(dimx);
        end
       
        stepsize = norm(t*desc_dir);
        dual_stepsize= norm(ynew-yk);
        xk = xnew;
        yk = ynew;
    end

    %generate output structure

    % compute complementary slack
    slack_cones =zeros(n_cones,1);
    for i=1:n_cones
        slack_cones(i) = abs(dot(gxk(block_begin(i):block_end(i)) ...
                                     ,yk(block_begin(i):block_end(i))));
    end
    comp_slack = max(slack_cones);
    % verify feasibility
    [gxk,~] = gj(xk); 
    gxk = gxk(:);    
    feasibility=max(abs(min(spectral_decomposition(gxk,mj), 0)));

    x= xk;
    [fval,grad] = fun(x); 
    output.iterations = k;
    output.cputime = toc;
    output.constrviolation = feasibility;
    output.firstorderopt = norm_lag;
    output.compslack = comp_slack;
    output.bestfeasible.x = x;
    output.bestfeasible.fval = fval;
    output.bestfeasible.constrviolation = output.constrviolation;
    output.bestfeasible.firstorderopt = output.firstorderopt;
    lambda.gj = yk;
    hessian = matB;
    if k==0
        output.stepsize=[];
    else
        output.stepsize = stepsize;
    end    
    if norm_lag <options.OptimalityTolerance
        % stop if the norm of the lagrangian is small
        exitflag = 1;
        output.message = append('Convergence by norm(Grad_Lag)<', ...
            num2str(options.OptimalityTolerance));
        disp(output.message);
    elseif ~is_linear_system_solved
        % stop if maximum iterations is reached
        exitflag = 0; 
        output.message = 'Linear System too stiff, algorithm cannot continue';
        disp(output.message);    
    elseif norm_da < options.StepTolerance
        % stop if norm of the direction vector is too small
        exitflag = 2;
        output.message = append('Convergence by norm of direction vector da <',...
            num2str(options.StepTolerance));
        disp(output.message);
    elseif stepsize < options.StepTolerance
        % stop if stepsize is too small
        exitflag = 2;
        output.message = append('Convergence by step size: stepsize <',...
            num2str(options.StepTolerance));
        disp(output.message);
    elseif k==options.MaxIterations 
        % stop if maximum iterations is reached
        exitflag = 0; 
        output.message = 'Number of iterations exceeded options.MaxIterations';
        disp(output.message);
    else
        warning('fdipa: unexpected algorithm exit.')
    end
end