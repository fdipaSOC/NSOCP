function optobj = fdipa_options(varargin)
% construct an object containing all the argument to call fdipa(..) 
% INPUT:
%   - varargin: Any number of arguments, separated by comma, indicating the 
%               options used in the execution of fdipa. the format is:
%               'nameOption1',valueOption1,'nameOption2',valueOption2,...
%               Valid options are the following: Display, SpecifyConstraintGradient, 
%               SpecifyObjectiveGradient, MaxIterations, ConstraintTolerance, 
%               OptimalityTolerance, StepTolerance, HessianFcn, HessianApproximation, 
%               ParXi, ParEta, ParNu, ParPhi, ParCI, ParCS, ParLambdam. See documentation for details.
% OUTPUT: 
%   - optobj: object of class optionsClass with the configuration for the execution of fdipa
%
% USAGE:
%   fdipa_options('MaxIterations',100,'StepTolerance',1e-12)
    
    % create an empty object of class options
    optobj = options_class();
    if nargin==0
        return
    end
    if mod(nargin,2)==1 
        error('fdipa_options:Wrong number of arguments');
    else
        for k=1:(nargin/2)
            switch varargin{2*(k-1)+1}
                case 'Display'
                    msj=optobj.edit('Display',varargin{2*k});
                case 'SpecifyConstraintGradient'   %ignored, gradient of constraint is required
                    msj=optobj.edit('SpecifyConstraintGradient',varargin{2*k});
                case 'SpecifyObjectiveGradient'      
                    msj=optobj.edit('SpecifyObjectiveGradient',varargin{2*k}) ;
                case 'MaxIterations'      
                    msj=optobj.edit('MaxIterations',varargin{2*k}) ;
                case 'ConstraintTolerance'       
                    msj=optobj.edit('ConstraintTolerance',varargin{2*k}) ;
                case 'OptimalityTolerance'       
                    msj=optobj.edit('OptimalityTolerance',varargin{2*k}) ;
                case 'StepTolerance'         
                    msj=optobj.edit('StepTolerance',varargin{2*k}) ;
                case 'HessianFcn'      
                    msj=optobj.edit('HessianFcn',varargin{2*k}) ;
                case 'HessianApproximation'      
                    msj=optobj.edit('HessianApproximation',varargin{2*k}) ;
                case 'ParXi'        
                    msj=optobj.edit('ParXi',varargin{2*k}) ;
                case 'ParEta' %Armijo parameter for the line search 
                    msj=optobj.edit('ParEta',varargin{2*k}) ;
                case 'ParNu'        
                    msj=optobj.edit('ParNu',varargin{2*k}) ;
                case 'ParPhi'       
                    msj=optobj.edit('ParPhi',varargin{2*k}) ;
                case 'ParCI'        
                    msj=optobj.edit('ParCI',varargin{2*k}) ;
                case 'ParCS'        
                    msj=optobj.edit('ParCS',varargin{2*k}) ;
                case 'ParLambdam'       
                    msj=optobj.edit('ParLambdam',varargin{2*k}) ;
                otherwise 
                    disp('The indicated variable does not belong to the set of options');
                    return      
            end
        end
        if msj==2 
            disp('fdipa_options:Wrong datatype');
            return
        elseif msj==3
            disp('fdipa_options:Parameter value out of Range');
            return
        end 
    end
 end
