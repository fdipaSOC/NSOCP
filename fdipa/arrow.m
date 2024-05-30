function arrz=arrow(z,mj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the arrow matrix of the vector. See [2] page 1324 for definition.
% When mj is given, returns a block diagonal matrix containing 
% arrow of (z(1), ...,z(mj(1)) ) in the first block and 
% arrow of (z(mj(1)+ ... + mj(k-1)+1), ..., z(mj(1)+ ... + mj(k))) 
% in the k-th block for 1< k < length(mj)
% 
% Computes the arrow function of the vector z = (z1,...,zJ) 
% of dimensions (m1,...,mj).
% INPUTS:
%  - z:  input vector in R x R^{m-1} for which we compute the arrow matrix.
%
% OPTIONAL INPUTS:
%  - mj: vector with the dimension of the cones when a block arrow matrix is needed.
%
% OUTPUTS:
%  - arrz: arrow matrix of the vector z, or block arrow matrix of z in the dimension mj.
%
% [2] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) 
% A feasible direction algorithm for nonlinear second-order cone programs, 
% Optimization Methods and Software, 34:6,1322-1341, DOI: 10.1080/10556788.2018.1506452
 
    % transpose z into a column vector if needed
    z=z(:);
    m=length(z);

    switch nargin
        case 1
            % construction of arrow matrix
            arrz = [z(1) z(2:m)';z(2:m) z(1)*eye(m-1)]; 
        case 2
            if sum(mj)~=m
                error('fdipa:arrow: Dimension missmatch,  mj(end) != len(z).');
            else
                % transpose mj into a column vector if needed
                mj = mj(:);
                n=length(mj);
                if n==1
                    arrz = arrow(z);
                else
                    arrz = zeros(m);
                    block_end = mj;
                    block_begin =ones(n,1);
                    for i=2:n
                        block_end(i)=block_end(i-1)+mj(i);
                        block_begin(i)=block_end(i-1)+1;
                    end
                    for i =1:n
                        arrz(block_begin(i):block_end(i),block_begin(i):block_end(i)) = arrow(z(block_begin(i):block_end(i)));
                    end
                end
            end
        otherwise
            error('fdipa:arrow: Wrong number of inputs'); 
    end
end