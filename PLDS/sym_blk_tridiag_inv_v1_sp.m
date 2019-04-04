function  [V, VV] = sym_blk_tridiag_inv_v1_sp(diagBlk, offDiagBlk, T) %scalar expression
%function  [D, OD] = sym_blk_tridiag_inv_v1_sp(AA,BB,adx,bdx)
% Compute block tridiagonal terms of the inverse of a *symmetric* block
% tridiagonal matrix. 
% 
% Note: Could be generalized to non-symmetric matrices, but it is not currently implemented. 
% Note: Could be generalized to compute *all* blocks of the inverse, but
% it's not necessary for my current application and hence not implemented.
%
% Note: AA and BB could be combined into a single variable, but I think that might get confusing. 
% 
% Input: 
%   AA - (n x n x Ka) unique diagonal blocks 
%   BB - (n x n x Kb) unique off-diagonal blocks
%  adx - (T x 1) vector such that (i,i)th block of A is    [1:T]'
%                   A_{ii} = AA(:,:,adx(ii))
%  bdx - (T-1 x 1) vector such that (i,i+1) block of sA is  ones(T-1,1)
%                   A_{i,i+1} = AA(:,:,bdx(ii))
% 
% Output: 
%   D  - (n x n x T) diagonal blocks of the inverse        smooth.V
%  OD  - (n x n x T-1) off-diagonal blocks of the inverse  smooth.VV
% 
% From: 
% Jain et al, 2006
% "Numerically Stable Algorithms for Inversion of Block Tridiagonal and Banded Matrices"
% (c) Evan Archer, 2014
%

%
% simplified version of sym_blk_tridiag_inv_v1, for a specific case of use.
% simplified version for scalar case Kalman smoother purpose
%
% change input parameter
% Input: 
%   diagBlk - (1 x T) unique diagonal blocks 
%   offDiagBlk - (scalar) unique off-diagonal blocks
%   T - (scalar) size of seq
% change output parameter
% Output: 
%   V  - (1 x T) diagonal blocks of the inverse        
%  VV  - (1 x T-1) off-diagonal blocks of the inverse    
% Sile Hu 2016-11-19
%

    offDiagBlk = -offDiagBlk; % we gotta make them the negative of the blocks --scalar

    S = zeros(T-1,1);
    V = zeros(T, 1); % diagonal  
    VV = zeros(T-1, 1); % off diagonal

    S(T-1) = offDiagBlk/ diagBlk(T);
    for idx = (T-2):-1:1
       S(idx) = offDiagBlk / (diagBlk(idx+1) - S(idx+1)*offDiagBlk');
    end
    % Compute diagonal and off-diagonal blocks
    V(1) = 1/(diagBlk(1) - offDiagBlk*S(1));
    VV(1) = S(1)*V(1);
    for idx = 2:T-1
       V(idx) = (diagBlk(idx) - offDiagBlk * S(idx))\(1 + offDiagBlk*V(idx-1)*S(idx-1)); 
       VV(idx) = S(idx)'*V(idx);
    end
    V(T) = (1 + offDiagBlk*V(T-1)*S(T-1))/diagBlk(T);
end