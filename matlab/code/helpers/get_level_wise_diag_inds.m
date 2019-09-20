function all_inds = get_level_wise_diag_inds(K, n_levels, varargin)
% all_inds = get_level_wise_diag_inds(K, n_levels)
%
% Find the linear indices that will give the diagonal elements of each
% 'level' in a tensor. The levels run along the last mode and all other
% modes must be diagonal, i.e. all modes in the tensor except the last
% mode must have the same size.
%
% K: size of all but the last mode (must be the same for all but the last
% mode).
%
% n_levels: size of the last mode.
%
% varargin{1}: number of modes on each level. Default: 2.

if isempty(varargin)
    n_modes = 2;
else
    n_modes = varargin{1};
end
    

diag_inds_first_level = get_diagonal_indices(K, n_modes);

inds_level_one = diag_inds_first_level; %(1:(K+1):K^2)';
to_add = 0:K^n_modes:((n_levels-1)*K^n_modes);
to_add = repelem(to_add, K);
all_inds = repmat(inds_level_one, [n_levels, 1]) + to_add';
end