function Xs_array = cell_array_to_nd_array(Xs_cell_array)
    n_tensors = length(Xs_cell_array);
    obs_dimensions = size(Xs_cell_array{1});
    n_observation_modes = length(obs_dimensions);
    
    
    cell_as_mat = cell2mat(Xs_cell_array);

    reshaped_mat = reshape(cell_as_mat, [obs_dimensions, n_tensors]);
    
    %permute_vector = [1:(n_observation_modes-1), n_observation_modes+1, n_observation_modes];
    permute_vector = [1:2, 4:(n_observation_modes+1), 3];

    permuted_mat_obs_in_last_mode = permute(reshaped_mat, permute_vector);
    
    Xs_array = permuted_mat_obs_in_last_mode;
end