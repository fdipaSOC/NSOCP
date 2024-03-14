function [lambda1,lambda2,u1,u2]=spectral_decomposition(w_vec,mj)
% Computes the spectral decomposition of the vector w_vec \in K^m 
% with with respect to the cone K^m, m = length(w_vec)). 
% For m >= 2. For m=1 it returns  a single vector u1, and with 
% corresponding lambda values lambda1 = lambda2 = w_vec(1).
%
% When mj is provided it computes the spectral decomposition on 
% different coordinates separately. In this case u1 and u2 are 
% the concatenation of all the spectral vectors and lambda1, lambda2_block
% are vectors with the spectral values on each cone.
%
% INPUTS:
%  - w_vec: vector in K^m1 x ...x K^mJ, where mj=(m1,... ,mJ), 
%        if the input mj is omitted the default is mj = length(w_vec).
% OPTIONAL INPUT: 
%  - mj: a vector with dimensions of the cones, if the spectral 
%        decomposition on each individual cone is needed.
%
% OUTPUTS:
%  - lambda1: first spectral value
%  - lambda2: second spectral value
%  - u1: first vector of the spectral decomposition
%  - u2: second vector of the spectral decomposition
%
% USAGE:
    if length(w_vec) == 1
        u1 = 1;
        u2 = 0;
        lambda1 = w_vec(1);
        lambda2 = w_vec(1);
        return
    end
    switch nargin
        case 1
            dim_m = length(w_vec);
            
            w1 = w_vec(1);
            w_bar = w_vec(2:end);
            if norm(w_bar)~=0
                % Spectral decomposition of w_vec
                v = 1/norm(w_bar)*w_bar;
                u1 = 0.5*[1;-v];
                u2 = 0.5*[1; v];
                lambda1 = w1 - norm(w_bar);
                lambda2 = w1 + norm(w_bar);
            
            else
                % if norm(w_bar) = 0 we can take any unitary v, so 
                % we choose v to be the first cannonical vector
                v = [1;zeros(dim_m-2,1)];
                u1 = 0.5*[1;-v];
                u2 = 0.5*[1; v];
                lambda1 = w1;
                lambda2 = w1;
            
            end
            return
        case 2
            block_begin = ones(length(mj),1);
            block_end = mj;
            for i=2:length(mj)
                block_end(i)=block_end(i-1)+mj(i);
                block_begin(i)=block_end(i-1)+1;
            end
                 
            lambda1 = zeros(length(mj),1);
            lambda2 = zeros(length(mj),1);
            u1 = zeros(length(w_vec),1);
            u2 = zeros(length(w_vec),1);
            
            for i= 1:length(mj)
                [lambda1_block,lambda2_block,u1_block,u2_block] = spectral_decomposition(w_vec(block_begin(i):block_end(i)));
                lambda1(i) = lambda1_block;
                lambda2(i) = lambda2_block;
                u1(block_begin(i):block_end(i))= u1_block;
                u2(block_begin(i):block_end(i))= u2_block;
            end 
        otherwise
            error('fdipa:spectral_decompositionrrow: Wrong number of inputs');
    end
end