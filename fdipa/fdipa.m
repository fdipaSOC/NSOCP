function [x,fval,exitflag,output,lambda,grad,hessian] = fdipa(varargin)
% FDIPA : feasible direction interior-point algorithm
% as presented in [1]. This is an extension to NSOCP of the feasible 
% direction interior-point algorithm (FDIPA) proposed by Herskovits 
% in [2] for smooth constrained nonlinear programming (NP).
%
% The algorithm follows closely [1] page 1330.
% [1] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) A feasible 
%     direction algorithm for nonlinear second-order cone programs, 
%     Optimization Methods and Software, 34:6, 1322-1341, 
%     DOI: 10.1080/10556788.2018.1506452
% [2] J. Herskovits, Feasible direction interior-point technique for 
%     nonlinear optimization, J. Optim. TheoryAppl. 99 (1998), pp. 121–146.
% 
% fdipa attempts to find a minimum of a scalar function with 
% second order cone constraints, starting at an initial estimate.
%
% INPUTS:
%  - fun:  Objective function, the gradient must be provided as a second output.
%  - x0:  Starting point for the optimization.
%  - [gj]: The funcion that computes the cone constraint gj. gj accepts a 
%         vector x and returns the vector gj as first output and grad gj as 
%         second out put. See documentation for details.
%
% OPTIONAL INPUTS:
%  - [mj]: If the problem has more than one cone constraint, then 
%         gj=[g1;...;gJ] and the dimension of each cone should be specified 
%         through the vector mj=[m1,...,mJ]. The gradients of the constraints 
%         must be supplied and the SpecifyConstraintGradient option must by 'on'. 
%  - y0:  An initial estimate of the lagrange multiplier assiciated with the 
%         cone constraints. This argument is optional, if you do not want to 
%         provide it set y0 = [], in this case the default value of y0 is 
%         (1, 0, 0, ... , 0) for each cone.
%  - options: Object of the class optionClass with the configurations for the 
%         execution of the algorithm.
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
%        problem.options: Options created with fdipaoptions
%
% OUTPUTS:
%  - x: argmin(fun)
%  - fval: Value of fun at the solution x.
%  - exitflag: Indicator of the status of the execution, has the value 0 
%              in successful.
%  - output: A structure with information about the optimization. The 
%            fields of the structure are the following:
%        output.iterations: Number of iterations taken
%        output.constrviolation: number of contraints not satisfied
%        output.stepsize: Length of last displacement in x
%        output.firstorderopt: Measure of first-order optimality
%        output.message: Exit message
%        output.bestfeasible: structure with the best feasible solution
%            found
%  - lambda: Structure with fields containing the Lagrange multipliers at the 
%        solution x:
%        lambda.gj:  lagrange multiplier corresponding to the inequality contraints
%  - grad: Gradient of fun at the solution x.
%  - hessian: Hessian approximation at the solution x.
%
% USAGE:
%   x = fdipa(fun,x0,[gj])
%   x = fdipa(fun,x0,[gj],[mj])
%   x = fdipa(fun,x0,[gj],[mj],y0)
%   x = fdipa(fun,x0,[gj],[mj],y0,options)
%   x = fdipa(problem)
%   
%   [x,fval] = fdipa(...)
%   [x,fval,exitflag] = fdipa(...)
%   [x,fval,exitflag,output] = fdipa(...)
%   [x,fval,exitflag,output,lambda] = fdipa(...)
%   [x,fval,exitflag,output,lambda,grad] = fdipa(...)
%   [x,fval,exitflag,output,lambda,grad,hessian] = fdipa(...)
%
% See ./documentation/doc_fdipa.pdf for details
% 

    t0 = cputime;
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
    [~,gradfx0] = fun(x0);   
    [gx0,gradgx0] = gj(x0); 
    dimx = length(x0);   
    if length(gradfx0)~=dimx
        error('fdipa: gradfx has incompatible dimensions.');
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
    
    % check if y0 is provided, if not we take the direction (1,0,...,0) on 
    % each cone, which is always feasible
    if min(size(y0))==0
        y0 = zeros(length(gx0),1);
        for i=1:n_cones
            y0(block_begin(i))=1;
        end
    end
    dimg = length(y0); 
      
    %check if the inputs have the correct dimensions
    if length(gx0)~=dimg
        error('fdipa: gj(x) has incompatible dimensions.');
    end 
    
    [nrows_gradg, ncol_gradg]=size(gradgx0);
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
        % case 'mod-newton'
        %     if ~isa(options.HessianFcn,'function_handle')
        %         error(append('fdipa: options.HessianApproximation = mod-newton, you ', ...
        %             'must supply an update of the Hessian or change ', ...
        %             'options.HessianApproximation'));
        %     else
        %         b_update = options.HessianFcn;
        %         matB=eye(dimx); 
        %     end
        % case 'on'
        %     if ~isa(options.HessianFcn,'function_handle')
        %         error(append('fdipa: options.HessianApproximation = on, you ', ...
        %             'must supply an update of the Hessian or change ', ...
        %             'options.HessianApproximation'));
        %     else
        %         b_update = options.HessianFcn;
        %         matB=eye(dimx); 
        %     end
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


    % cast vectors as column vectors    
    gx0 = gx0(:);
    
    % Start of the algorithm    
    % Step 0
    xk = x0; 
    yk = y0;
    matB=eye(dimx); 
    lamb_min =spectral_decomposition(gx0,mj);
    total_time = cputime-t0;    
    % Check feasible point, min(lamb_min) < 0 means that the point is unfeasible
    %if min(lamb_min)<-options.ConstraintTolerance(1)
    if min(lamb_min)<0
        x = x0;
        [fval,grad] = fun(x); 
        exitflag = -1;
        % output
        output.iterations = 0;
        output.cputime = total_time;
        output.constrviolation = sum(lamb_min<0);
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
        time_k = cputime;
        % Step 1
        % evaluate Lagrangian at (x,yk)
        [gxk,grad_gxk] = gj(xk); 
        gxk = gxk(:);       
        [fxk,grad_fxk] = fun(xk);
        grad_fxk = grad_fxk(:);
        norm_Lag = norm(grad_fxk -grad_gxk'*yk);  

        %Report iteration information
        % Start screen information table in case 'iter'             
        if strcmp(options.Display,'iter')
            % compute complementary slack
            gjyj =0;
            for i=1:n_cones
                gjyj = gjyj + abs(dot(gxk(block_begin(i):block_end(i)) ...
                                             ,yk(block_begin(i):block_end(i))));
            end
            % verify feasibility
            lamb_min=spectral_decomposition(gxk,mj);        
            if k==0
                fprintf(append(repelem(' ',32),'Feasibility  First-order    Dual step  Comp. Slack.\n'));
                fprintf(append('Iter         f(x)    Step-size    (lambmin)   optimality  ',...
                    '  |yk+1-yk|      |gj.yj|\n'));
                fprintf(['%-4d  %11.6e  ',repelem(' ',13),'%11.3e  %11.3e  ',repelem(' ',13),'%11.3e\n'],...
                    k,fun(xk), min(lamb_min), norm_Lag, gjyj);
            else
                fprintf('%-4d  %11.6e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e\n', ...  
                    k,fun(xnew),stepsize,min(lamb_min_t),...
                    norm_Lag,dual_stepsize, gjyj);
            end
        end
        if k == options.MaxIterations
            break
        end
   
        % stop if first order optimallity is satified
        if norm_Lag <options.OptimalityTolerance
            break 
        end
        if k>0 && stepsize < options.StepTolerance
            break
        end

        % (i) Solving 1st System 
        arrw_y=arrow(yk,mj); % Arrow matrix of yk
        arrw_g=arrow(gxk,mj); % Arrow matrix of gj(xk)        
        % solve for both (da, ya) simultaneously and split it later
        daya=[matB -grad_gxk'; arrw_y*grad_gxk arrw_g]\[-grad_fxk; zeros(sum(mj),1)];
        da = daya(1:dimx);
        ya = daya((dimx+1):(dimx + sum(mj)));
        norm_da = norm(da); 

        %  stop if norm of the direction vector is too small
        if norm_da < options.StepTolerance
            break
        end 

        % (ii) Solving 2nd system 
        
        % solve for both (db, yb) simultaneously and split it later
        dbyb = [matB,-grad_gxk';arrw_y*grad_gxk, arrw_g]\[zeros(dimx,1);yk];
        db = dbyb(1:dimx); 
        yb = dbyb((dimx+1):end);

        % (iii) Compute rho
        
        pi_dGf = db'*grad_fxk;
        if pi_dGf>0
            rho = min(options.ParPhi*norm_da^2,...
                 (options.ParXi-1)*(da'*grad_fxk)/pi_dGf); 
        else
            rho = options.ParPhi*norm_da^2; 
        end

        % (iv) Compute feasible descent direction
        desc_dir = da + rho*db;  
        yk_hat = ya + rho*yb;
        
        
        % Step 2: Armijo linear search
        t=1;
        xtd=xk + t*desc_dir;
        fxt=fun(xtd);
        gxt=gj(xtd);
        lamb_min_t=spectral_decomposition(gxt,mj);
        iter_lin=0;
        
        % to accept a step size we check if with the choosed step we have 
        % enough descend, the point we reach is feasible and an additional 
        % technical spectral condition specified in the algorithm
        has_enough_descend = fxt-fxk<=t*options.ParEta*desc_dir'*grad_fxk;
        %we require that we are strictly in the interior since it is a
        %interior point algorithm
        is_feasible = min(lamb_min_t)>0;
        %is_feasible = min(lamb_min_t)>-options.ConstraintTolerance(1);
        [lambda1_gxk,lambda2_gxk,u1_gxk,u2_gxk] ...
            = spectral_decomposition(gxk,mj);
        % spectral conditions
        satisfy_spectral_condition = true;
        [lambda1_gxt,lambda2_gxt] ...
            = spectral_decomposition(gxt,mj);
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
        
        while (~has_enough_descend)||(~is_feasible) ...
                ||(~satisfy_spectral_condition)
            t=t*options.ParNu; 
            fxt=fun(xk + t*desc_dir);
            gxt=gj(xk + t*desc_dir); 
            lamb_min_t=spectral_decomposition(gxt,mj);
            iter_lin=iter_lin+1;
            has_enough_descend = fxt<=fxk+t*options.ParEta*desc_dir'*grad_fxk;
            %is_feasible = min(lamb_min_t)>-options.ConstraintTolerance(1);
            is_feasible = min(lamb_min_t)>0;
            % spectral conditions
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
                % update y
                ynew(block_begin(i))= max(ya(block_begin(i)),options.ParCI);
            else
                % Projection of yak in the subspace u1,u2
                % multiply by 2 because u1,u2 are not unit vectors
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
        %check for y feasible?
        if min(spectral_decomposition(ynew,mj))<0
            error("multiplier outside the cone")
        end

        
        % Update Bk
        if mod(k+1,options.HessianResetIterations)==0
            matB = eye(dimx);
            %fprintf('reset iteration')
        else
            matB = b_update(xnew,xk,ynew,yk,fun,gj, matB);
        end 
        % A quantifiable way of checking that Assumtion 3.5 in [1] is satisfied 
        % is to reset the matrix B to the identity if the condition number is too small
        % Condition on sigma2 not being too large
        if eigs(matB,1) > options.ParSigma2 
            matB = eye(dimx);
            fprintf('Sigma2 too large\n')
        end
        % Condition on the ratio sigma2/sigma1 not being too large
        if cond(matB) > options.NumericalConditioning %|| rcond(matB) < 1/options.NumericalConditioning
            matB = eye(dimx);
        end

        
        stepsize = norm(t*desc_dir);
        dual_stepsize= norm(ynew-yk);
        xk = xnew;
        yk = ynew;

        delta_time_k = cputime-time_k;
        total_time = total_time + delta_time_k;
    end

    %generate output structure
    delta_time_k = cputime-time_k;
    total_time = total_time+delta_time_k;
    x= xk;
    [fval,grad] = fun(x); 
    output.iterations = k;
    output.cputime = total_time;
    output.constrviolation = sum(lamb_min<0);
    output.firstorderopt = norm_Lag;
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
    if norm_Lag <options.OptimalityTolerance
        % stop if the norm of the lagrangian is small
        exitflag = 1;
        output.message = append('Convergence by norm(Grad_Lag)<', ...
            num2str(options.OptimalityTolerance));
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