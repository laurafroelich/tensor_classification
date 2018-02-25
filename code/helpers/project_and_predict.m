function predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us)

x_train_projected = tensor_projection(x_train, Us);
x_test_projected = tensor_projection(x_test, Us);

n_train_obs = length(y_train);
n_test_obs = length(y_test);

x_train_proj_mat = reshape(cell2mat(x_train_projected)', [ncomps^2, n_train_obs])';
x_test_proj_mat = reshape(cell2mat(x_test_projected)', [ncomps^2, n_test_obs])';

mnr_mod = mnrfit(x_train_proj_mat, y_train);
predictions = mnrval(mnr_mod, x_test_proj_mat);
end