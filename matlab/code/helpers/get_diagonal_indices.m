function diagonal_indices = get_diagonal_indices(k)
% diagonal_indices = get_diagonal_indices(k)
%
% Get linear indices of elements on multi-dimensional diagonal.
%
% Input: 
%      k : integer, dimensionality of multidimensional cube.
%
% Output:
%      vector containing linear indices of diagonal

        diagonal_indices = diag(cumsum([1:(k+1):k^2; k^2.*ones(k-1,k)]));
end