function Xs_array = cell_array_to_nd_array(Xs_cell_array)
% Xs_array = cell_array_to_nd_array(Xs_cell_array)
%
% Arrange n-mode arrays from cell array in (n+1)-mode array.
%
% Each n-mode array in the input 'Xs_cell_array' must have the same
% dimensions. An n-mode array with observations along the (n+1)'th mode
% is returned.
%
% Input:
%
% Xs_cell_array: Cell array containing numeric arrays of equal dimensions.
%
% Output:
%
% Xs_array: Numeric array with observations along last dimension.
    n_tensors = length(Xs_cell_array);
    obs_dimensions = size(Xs_cell_array{1});
    n_observation_modes = length(obs_dimensions);
    
    
    cell_as_mat = cell2mat(Xs_cell_array);

    reshaped_mat = reshape(cell_as_mat, [obs_dimensions, n_tensors]);
    
    permute_vector = [1:2, 4:(n_observation_modes+1), 3];

    permuted_mat_obs_in_last_mode = permute(reshaped_mat, permute_vector);
    
    Xs_array = permuted_mat_obs_in_last_mode;
end
