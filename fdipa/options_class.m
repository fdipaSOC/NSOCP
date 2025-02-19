classdef options_class < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CLASS OPTIONS
    properties 
        Display; %indicates how much information will be shown
        MaxIterations; %maximum number of iteration for the algorithm
        %ConstraintTolerance; %tolerance in the values of the constrain function
        OptimalityTolerance; % stopping test: Lagragian
        StepTolerance; % stopping test: da - direction 
        LowerOptimalityBound; % stopping test, value of objective function
        HessianFcn; % a function handle to custom Hessian function
        HessianApproximation; % indicates what hessian approximation will be used
        ParXi;
        ParEta;
        ParNu;
        ParPhi;
        ParCI;
        ParCS;
        ParLambdam;
        ParSigma1;
        ParSigma2;
        NumericalConditioning;
        LinearSystemTolerance;
    end
   
    methods
        function opt = options_class() 
        %Constructor of the class, set the default values of the parameters
            opt.Display = 'final';
            opt.MaxIterations = 1000;
            %opt.ConstraintTolerance = 1e-6;
            opt.OptimalityTolerance = 1e-6;
            opt.StepTolerance = 1e-10;
            opt.LowerOptimalityBound = -realmax; % stopping test, value of objective function
            opt.HessianFcn  = NaN;
            opt.HessianApproximation= 'default';
            opt.ParXi= 0.7;
            opt.ParEta= 1e-4;
            opt.ParNu= 0.7;
            opt.ParPhi= 1;
            opt.ParCI= 1e-9;
            opt.ParCS = 1e9;
            opt.ParLambdam = 1e-4;
            opt.ParSigma1 = 1e-15;
            opt.ParSigma2 = 1e15;
            opt.NumericalConditioning = 1e15;
            opt.LinearSystemTolerance = 1e-4;
        end
      
        function status = edit(opt,varName,value) 
        % Method that changes the value of one property of the class at a time
        % 
        % INPUTS:
        %    opt: option_class object, used for the declaration of the method
        %    varName: name of the parameter whose value needs to be changed
        %    value: value of the parameter we want to modify.
        %
        % OUTPUTS:
        %    status: Returns 0 if the chosen variable do not belong to the set of admisible options
        %            Returns 1 if the variable belongs to the set of admisible options and if 
        %            chosen value is acceptable it changes the value of such parameter.
        %            Returns 2 if the datatype of 'value' is incorrect. 
        %            Returns 3 if the variable belongs to the set of admisible options and but 
        %            chosen value is is not admisible for such parameter
        %
            status=1;
            switch varName
                case 'Display'
                    if ischar(value)
                        switch value
                            case 'off' 
                                opt.Display = value;
                            case 'none' 
                                opt.Display = value;
                            case 'iter' 
                                opt.Display = value;
                            case 'final' 
                                opt.Display = value;
                            otherwise
                                status=3;
                        end
                    else
                        status=2;
                    end
                case 'MaxIterations'
                    if isnumeric(value)
                        if value>0
                            opt.MaxIterations = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                % case 'ConstraintTolerance'
                %     if isnumeric(value)==1
                %         if value>0
                %             opt.ConstraintTolerance = value;
                %         else
                %             status=3;
                %         end
                %     else
                %         status=2;
                %     end
                case 'OptimalityTolerance'
                    if isnumeric(value)
                        if value>0
                            opt.OptimalityTolerance = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'StepTolerance'
                    if isnumeric(value)
                        if value>0
                            opt.StepTolerance = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end 
                case 'LowerOptimalityBound'
                        if isnumeric(value)
                            opt.LowerOptimalityBound = value;
                        else
                            status=2;
                        end 
                case 'HessianFcn'
                    if isa(value,'function_handle')
                        opt.HessianFcn = value;
                    else
                        status=2;
                    end
                case 'HessianApproximation'
                    if ischar(value)
                        switch value
                            case 'off' 
                                opt.HessianApproximation = value;
                            case 'default' 
                                opt.HessianApproximation = value;
                            case 'user-supplied' 
                                opt.HessianApproximation = value;
                            case 'bfgs' 
                                opt.HessianApproximation = value;
                            otherwise
                                status=3;
                        end
                    else
                        status=2;
                    end                    
                case 'ParXi'
                    if isnumeric(value)
                        if (value>0)&&(value<1)
                            opt.ParXi = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'ParEta'
                    if isnumeric(value)
                        if (value>0)&&(value<1)
                            opt.ParEta = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'ParNu'
                    if isnumeric(value)
                        if (value>0)&&(value<1)
                            opt.ParNu = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'ParPhi'
                    if isnumeric(value)
                        if value>0
                            opt.ParPhi = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'ParCI'
                    if isnumeric(value)
                        if value>0
                            opt.ParCI = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'ParCS'
                    if isnumeric(value)
                        if value>0
                            opt.ParCS = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'ParLambdam'
                    if isnumeric(value)
                        if min(value)>0
                            opt.ParLambdam = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'ParSigma1'
                        if isnumeric(value)
                            if min(value)>0
                                opt.ParSigma1 = value;
                            else
                                status=3;
                            end
                        else
                            status=2;
                        end 
                case 'ParSigma2'
                    if isnumeric(value)
                        if min(value)>0
                            opt.ParSigma2 = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end 
                case 'NumericalConditioning'
                    if isnumeric(value)
                        if min(value)>0
                            opt.NumericalConditioning = value;
                        else
                            status=3;
                        end
                    else
                        status=2;
                    end
                case 'LinearSystemTolerance'
                        if isnumeric(value)
                            if min(value)>0
                                opt.LinearSystemTolerance = value;
                            else
                                status=3;
                            end
                        else
                            status=2;
                        end  
                otherwise
                    status=0;
            end
        end
   end 
end 

