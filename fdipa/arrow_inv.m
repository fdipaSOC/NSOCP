function arrz=arrow_inv(z,mj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package has been downloaded from https://github.com/fdipaSOC/NSOCP
% This example is included in the as an application of the algorithm described 
% in [1]. See README.md for details.
% [1] Alfredo Canelas, Miguel Carrasco, Julio Lopez, Esteban Paduro (2024)
%     FDIPA-SOC: A MATLAB Package for Nonlinear Second-Order Cone Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the inverse of the arrow matrix of the vector. See [2] page 1324 for definition.
% When mj is given, returns a block diagonal matrix containing 
% inverse arrow matrix of (z(1), ...,z(mj(1)) ) in the first block and 
% inverse arrow matrix of (z(mj(1)+ ... + mj(k-1)+1), ..., z(mj(1)+ ... + mj(k))) 
% in the k-th block for 1< k < length(mj)
% 
% Computes the arrow function of the vector z = (z1,...,zJ) 
% of dimensions (m1,...,mj).
% INPUTS:
%  - z:  input vector in R x R^{m-1} for which we compute the inverse of the arrow matrix.
%
% OPTIONAL INPUTS:
%  - mj: indices of the dimension when a block arrow matrix is needed.
%
% OUTPUTS:
%  - arrz: arrow matrix of the vector z, or block arrow matrix of z in the dimension mj.
%
% [2] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) 
% A feasible direction algorithm for nonlinear second-order cone programs, 
% Optimization Methods and Software, 34:6,1322-1341, DOI: 10.1080/10556788.2018.1506452
    z=z(:);
    m=length(z);
    switch nargin
        case 1
            % construction of the inverse of the arrow matrix
            det_z=z(1)^2-norm(z(2:m))^2;
            arrz(1,1)=z(1)/det_z;
            arrz(1,2:m)=-z(2:m)'/det_z;
            arrz(2:m,1)=-z(2:m)/det_z;
            arrz(2:m,2:m)=(eye(m-1)+z(2:m)*z(2:m)'/det_z)/z(1);
            
        case 2
            if sum(mj)~=m
                error('fdipa:arrow_inv: Dimension missmatch.');
            else
                % convert mj in a column vector if needed 
                mj=mj(:);
                n=length(mj);
                if n==1
                    arrz = arrow_inv(z);
                else
                    arrz = zeros(m);
                    block_end = mj;
                    block_begin =ones(n,1);
                    for i=2:n
                        block_end(i)=block_end(i-1)+mj(i);
                        block_begin(i)=block_end(i-1)+1;
                    end
                    % the inverse of the arrow matrix is computed block-wise
                    for i =1:n
                        arrz(block_begin(i):block_end(i),block_begin(i):block_end(i)) = arrow_inv(z(block_begin(i):block_end(i)));
                    end
                end
            end
        otherwise
            error('fdipa:arrow_inv: Wrong number of inputs'); 
    end
end
