function diagonal_indices = get_diagonal_indices(k, varargin)
% diagonal_indices = get_diagonal_indices(k)
%
% Get linear indices of elements on multi-dimensional diagonal.
%
% Input: 
%      k : integer, dimensionality of multidimensional cube.
%
% varargin{1}: number of modes on each level. Default: 2.
%
% Output:
%      vector containing linear indices of diagonal

if isempty(varargin)
    n_modes = 2;
else
    n_modes = varargin{1};
end

if n_modes == k
    diagonal_indices = diag(cumsum([1:(k+1):k^2; k^2.*ones(k-1,k)]));
end

if n_modes == 2
    diagonal_indices = (1:(k+1):k^2)';
end

end