function data_cell = mat_to_cell(data_mat)
% assumes observations are along last dimension

size_data = size(data_mat);
n_dims = length(size_data);
n_obs = size_data(end);
data_cell = cell(n_obs, 1);
data_mat_perm = permute(data_mat, [n_dims, 1:(n_dims-1)]);

for iobs = 1:n_obs
    data_cell{iobs} = reshape(data_mat_perm(iobs,:), size_data(1:end-1));
end
end