%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example of fdipa
% Example 2: Quadratic objective with rotated conic constraint
%
% $$\min f(x,y,s) =  s_1^2 + \cdots + s_m^2+ x^2 + y^2 - 2 s_1 s_2 + x + y +2 s_1+s_2,$$
% Subject to
% $$\|s\|^2 \leq 2 x y, ~x\geq 0, ~y\geq 0,$$
%
% The constraint can be writen as a second order cone constraint using
%
% $$
% \left\| \left[ \begin{array}{c} \sqrt{2}s\\ x-y\end{array}\right] \right\| \leq x+y \Leftrightarrow 
% \left(\begin{array}{ccc}
% 1 & 1 & 0\\
% 0 & 0 & \sqrt{2}I_m \\
% 1 & -1 & 0 \\
% \end{array}\right)
% \left(\begin{array}{c} x\\y \\ s \end{array}\right) \in \mathcal{K}^{m+2}.
% $$
%
%% Input Functions
% 
% First we need to define the objective function and its gradient 
% as a MATLAB(R) function, i.e., here m is a parameter that indicates the 
% size of the vector s
%
%    function [fun,grad_f]=fun_rotated_conic(x,m)
%    % Quadratic objective function for the example with rotated conic constraint
%        x= x(:);
%        fun= sum(x.^2)-2*x(3)* x(4) + x(1)+x(2)+2*x(3)+x(4);
%        grad_f=2*x+[1;1;2-2*x(4);1-2*x(3);zeros(m-2,1)] ;
%    end
% 
% In equivalent way to the restrictions.
%
%    function [g_fun,grad_g]=g_rotated_conic(x,m)
%    % contraint function for the rotated conic constraint
%        x=x(:);
%        matA = zeros(m+2,m+2);
%        matA(1,1) = 1;
%        matA(1,2) = 1;
%        matA(m+2,1) = 1;
%        matA(m+2,2) = -1;
%        matA(2:(m+1),3:m+2) = eye(m);
%        g_fun = matA * x;
%        grad_g = matA;
%    end                           
%
%% Options
% 
% To create and modify the options of the algorithm, we only want to see a 
% final summary of the algorithm, this is done by setting |Display|  
% to |final| and all other settings as default, which is done using
%
%   myOptions = fdipa_options('Display','final');
%
%% Other Inputs Arguments
%
% Since the constraint can be written as a single conical constraint, 
% we can leave the parameter |mj = []|. 
% Finally you need to select a starting point |x0| which is done using 
% a construction that guarantees that the starting point is strictly feasible.
% Lastly, for the Lagrange multiplier |y0| we the default value.



seed = RandStream('mt19937ar','Seed',1);
%m=10;
%m=30;
%m=100;
m=1000;
%m=10000;

% we construct a feasible starting point using the geometry of the problem
a = 10*rand(seed);
b = 10*rand(seed);
c = -1+ 2*rand(seed,m,1);
x0  = [a;b; rand(seed)*sqrt(2 *a*b)/norm(c)*c];
% or we can use the subroutine to look for a feasible point using fdipa
%x0 = searchStartingPoint(m+2,@(x)g_rotated_conic(x,m),[]);


my_options = fdipa_options('Display','final','ParEta',0.5);
[~,~,~,output] = fdipa(@(x)fun_rotated_conic(x,m),x0,@(x)g_rotated_conic(x,m),...
    [],[],my_options);

% output for the paper 
fprintf('%d & %d & %11.5e & %3.2f \n',m, output.iterations, output.firstorderopt, output.walltime)

clear 'seed' 'a' 'b' 'c' 'x0' 'myoptions' 'fval' 'output' 'm' 'my_options'

function [fun,grad_f]=fun_rotated_conic(x,m)
% Quadratic objective function for the example with rotated conic constraint
    x= x(:);
    fun= sum(x.^2)-2*x(3)* x(4) + x(1)+x(2)+2*x(3)+x(4);
    grad_f=2*x+[1;1;2-2*x(4);1-2*x(3);zeros(m-2,1)] ;
end

function [g_fun,grad_g]=g_rotated_conic(x,m)
% contraint function for the rotated conic constraint
    x=x(:);
    matA = zeros(m+2,m+2);
    matA(1,1) = 1;
    matA(1,2) = 1;
    matA(m+2,1) = 1;
    matA(m+2,2) = -1;
    matA(2:(m+1),3:m+2) = sqrt(2)* eye(m);
    g_fun = matA * x;
    grad_g = matA;
end