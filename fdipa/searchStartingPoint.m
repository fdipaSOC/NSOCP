function [x0,exitflag] = searchStartingPoint(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutine to look for a initial feasible point for a second order cone program
% utilizing FDIPA-SOC
% INPUTS:
%  - dimx: Required input. Denotes a positive integer indicating the dimension 
%         n of the vector x in Problem.
%  - gj:  The function that computes the cone constraint gj. gj accepts a 
%         vector x and returns a vector gj with the evaluation of the constraint 
%         function as first output and the gradient of the constraints as 
%         second output. See documentation for details.
%
% OPTIONAL INPUTS:
%  - mj:  If the problem has more than one cone constraint, then 
%         gj=[g1;...;gJ] and the dimension of each cone should be specified 
%         through the vector mj=[m1,...,mJ].  
%  - xguess:  An initial guess to start the search of a feasible point.
% OUTPUTS:
%  - x0: argmin(fun)
%  - exitflag: Indicates if the x0 satisfy the constraints. 1 if it does, 0 if it does not.

    %parse an check the inputs
    if ~any([2,3,4]==nargin)
        error('fdipa: Wrong number of inputs.');
    end 
    dimx = varargin{1};
    gj = varargin{2};

    if nargin==4
        xguess = varargin{4};
    else
        xguess = [];
    end
    
    if min(size(xguess))==0
        xguess = zeros(dimx,1);
    end

    [gx,~] = gj(xguess);

    if nargin>2
        mj = varargin{3};
    else
        mj = [];
    end

    if min(size(mj))==0
        mj = length(gx);
    end

    %if the guess is feasible, we return it.
    lambda1 = spectral_decomposition(gx,mj);
    if min(lambda1)>0
        x0 = xguess;
        exitflag =  (min(lambda1)>0);
        return
    end
    %we add a small value since we need a strictly feasible point
    t_candidate = max(-lambda1) + 1e-2;
    x0_aux = [xguess;t_candidate];

    myoptions_aux = fdipa_options('Display','off','LowerOptimalityBound', 0);
    x_aux = fdipa(@(z) f_aux(z),x0_aux,@(z) g_aux(z,gj,mj),mj,[],myoptions_aux);
    x0 = x_aux(1:(end-1));
    [gx0,~] = gj(x0);
    lamb_min = spectral_decomposition(gx0,mj);
    exitflag =  (min(lamb_min)>0);
    if ~exitflag
        warning('searchStartingPoint: search of feasible starting point failed.')
    end
end

function [fx,grad_fx]=f_aux(x)
    fx = x(end);
    grad_fx = [zeros(length(x)-1,1);1];
end

function [gx,grad_gx]=g_aux(x,gj,mj)
    [gx,grad_gx] = gj(x(1:(end-1)));

    n_cones = length(mj); % number of cones
    block_begin = ones(n_cones,1); % index of the first coordinate of i-th cone
    block_end = mj; % index of the last coordinate of i-th cone 
    if n_cones>1
        for i=2:n_cones
            block_end(i)=block_end(i-1)+mj(i);
            block_begin(i)=block_end(i-1)+1;
        end
    end 

    grad_gx= [grad_gx,zeros(sum(mj),1)];

    for k=1:n_cones
        gx(block_begin(k)) = gx(block_begin(k)) + x(end);
        grad_gx(block_begin(k),:) = grad_gx(block_begin(k),:)+ [zeros(1,length(x)-1), 1];
    end
end