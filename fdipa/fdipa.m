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
%        output.constrviolation: Maximun of constraint functions
%        output.stepsize: Length of last displacement in x
%        output.firstorderopt: Measure of first-order optimality
%        output.message: Exit message
%        output.bestfeasible: structure with the best feasible solution
%            found
%  - lambda: Structure with fields containing the Lagrange multipliers at the 
%        solution x
%  - grad: Gradient of fun at the solution x.
%  - hessian: Hessian of fun at the solution x.
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
    % verification for the number of arguments is correct
    exitflag = 0;

    if ~any([1,3,4,5,6]==nargin)
        error('fdipa: Wrong number of inputs.');
    end


    if nargin==1
        if ~isstruct(varargin{1})
            error('fdipa: The input has not the correct format');
        end
        problem = varargin{1};
        if ~all(isfield(problem,{'objective','x0','gj'}))
            error(['fdipa: The input has not the correct format, ' ...
                'missing fields']);
        end
        fun = problem.objective;
        x0 = problem.x0;
        gj = problem.gj;
        if isfield(problem,'mj')
            mj = problem.mj;
        else
            mj = []
        end
        if isfield(problem,'mj')
            y0 = problem.y0;
        else
            y0 = []
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
  
    
    % check if the gradient of the objective function is provided, if not, stop.
    if nargout(fun)==1
        error('fdipa: You must supply the gradient of the objetive function.');
    end
    
    % check if the gradient of the constraints gj provided, if not, stop
    if nargout(gj)==1
        error('fdipa: You must supply the gradient of the cone constraints.');
    end
    
    % set the initial condition to a column vector
    x0=x0(:);
    
    % check compatibility of the input variables and functions
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
    switch options.HessianApproximation
        case 'off'
            matB=eye(dimx);
        case 'bfgs'
            matB=eye(dimx);
            b_update = @bfgs_update;
        case 'mod-newton'
            if ~isa(options.HessianFcn,'function_handle')
                error(append('fdipa: options.HessianApproximation = mod-newton, you ', ...
                    'must supply an update of the Hessian or change ', ...
                    'options.HessianApproximation'));
            else
                b_update = options.HessianFcn;
                matB=eye(dimx); 
            end
        case 'on'
            if ~isa(options.HessianFcn,'function_handle')
                error(append('fdipa: options.HessianApproximation = on, you ', ...
                    'must supply an update of the Hessian or change ', ...
                    'options.HessianApproximation'));
            else
                b_update = options.HessianFcn;
                matB=eye(dimx); 
            end
        case 'user-supplied'
            if ~isa(options.HessianFcn,'function_handle')
                error(append('fdipa: options.Hessian = supplied, you ', ...
                    'must supply an update of the Hessian or change ', ...
                    'options.HessianApproximation'));
            else
                b_update = options.HessianFcn;
                matB=eye(dimx); 
            end
        otherwise 
            error('fdipa:Choosen option for the Hessian is not valid.');
    end

    % cast vectors as column vectors    
    gx0 = gx0(:);
    
    % Start of the algorithm    
    % Step 0
    xk = x0; 
    yk = y0;
 
    lamb_min =spectral_decomposition(gx0,mj);
    total_time = cputime-t0;    
    % Check feasible point, min(lamb_min) < 0 means that the point is unfeasible
    if min(lamb_min)<-options.ConstraintTolerance(1)
        x = x0;
        [fval,grad] = fun(x); 
        exitflag = -1;
        % output
        output.iterations = 0;
        output.cputime = total_time;
        output.constrviolation = sum(lamb_min<=-options.ConstraintTolerance(1));
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
        % constraint and gradient evaluation
        [gxk,grad_gxk] = gj(xk); 
        % cast gxk as column vector
        gxk = gxk(:);
        % compute the first spectral value on each cone
        lamb_min=spectral_decomposition(gxk,mj);
            arrw_inv_g=arrow_inv(gxk,mj); % Inverse arrow matrix of gj(xk)
        arrw_y=arrow(yk,mj); % Arrow matrix of yk
        arrw_g=arrow(gxk,mj); % Arrow matrix of gj(xk)
            [fxk,grad_fxk] = fun(xk); % objective and gradient function evaluation
        grad_fxk = grad_fxk(:);
        
        % (i) Solving 1st System 
        % here we use that 'ya' can be solved explicitely in terms of 'da'
        matM1=arrw_inv_g*arrw_y*grad_gxk;
        da=-(matB+grad_gxk'*matM1)\grad_fxk;
        ya=-matM1*da;
        norm_da = norm(da); 
        % gradient of Lagrangian in (x,ya)
        grad_Lag_x=grad_fxk -grad_gxk'*ya; 
        % norm gradient of Lagrangian in (x,ya)
        norm_Lag = norm(grad_Lag_x);  
      
        %Report iteration information
        % Start screen information table in case 'iter'             
        if strcmp(options.Display,'iter')
            % compute complementary slack
            gjyj =0;
            for i=1:n_cones
                gjyj = gjyj + abs(dot(gxk(block_begin(i):block_end(i)) ...
                                             ,yk(block_begin(i):block_end(i))));
            end        
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
        %  stop if norm of the direction vector is too small
        if norm_da < options.StepTolerance
            break
        end 
        % stop if first order optimallity is satified
        if norm_Lag <options.OptimalityTolerance
            break 
        end
        if k>0 && stepsize < options.StepTolerance
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
        is_feasible = min(lamb_min_t)>-options.ConstraintTolerance(1);
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
            is_feasible = min(lamb_min_t)>-options.ConstraintTolerance(1);
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
        
        % New dual point
        % Spectral decomposition of g(xnew)
        [~,grad_gxk1] = gj(xnew);
        %[gxnew,grad_gxk1] = gj(xnew);
        %gxnew = gxnew(:);
        
        % update y step
        ynew = zeros(dimg,1);
        [~,~,u1_gxt,u2_gxt] = spectral_decomposition(gxt,mj);
        
        for i=1:n_cones
            if mj(i)==1
                % update y
                ynew(block_begin(i))= max(ya(block_begin(i)),0.001*norm(da)^2);
            else
                % Projection of yak in the subspace u1,u2
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
        grad_Lag_x=grad_fxk -grad_gxk'*ynew;
        [~,grad_fxk1] = fun(xnew);
        grad_fxk1 = grad_fxk1(:);
        grad_xnew=grad_fxk1 -grad_gxk1'*ynew;
        if ~strcmp(options.HessianApproximation,'off')
            if strcmp(options.HessianApproximation, 'mod-newton') 
                xkyk = zeros(dimx+dimg,1);
                xkyk(1:dimx) = xk;
                xkyk(dimx+(1:dimg)) = yk;
                matB = b_update(xkyk);
            elseif strcmp(options.HessianApproximation, 'bfgs') 
                %reset the Hessian after dimx iterations to ensure
                %Assumption 3.5 is satisfied
                if mod(k+1,dimx)==0
                    matB = eye(dimx);
                else
                    matB = b_update(xnew,xk,grad_xnew, grad_Lag_x, matB);
                end             
            else
                matB = b_update(xnew,xk,grad_xnew, grad_fxk, matB);
            end
            %%why is this here?
            %if sum(isnan(matB)>0)
            %    matB = eye(dimx);
            %end
        else
            matB = eye(dimx);
        end
        
        stepsize = norm(t*desc_dir);
        dual_stepsize= norm(ynew-yk);
        xk = xnew;
        yk = ynew;

        delta_time_k = cputime-time_k;
        total_time = total_time + delta_time_k;
    end

    %generate output
    delta_time_k = cputime-time_k;
    total_time = total_time+delta_time_k;
    x= xk;
    [fval,grad] = fun(x); 
    output.iterations = k;
    output.cputime = total_time;
    output.constrviolation = sum(lamb_min<-options.ConstraintTolerance(1));
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