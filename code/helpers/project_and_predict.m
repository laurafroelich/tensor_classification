function predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us, varargin)

if isempty(varargin)
    parafac_structure = false;
else
    parafac_structure = varargin{1};
end


x_train_projected = tensor_projection(x_train, Us);
x_test_projected = tensor_projection(x_test, Us);

n_train_obs = length(y_train);
n_test_obs = length(y_test);

if ~parafac_structure
    if iscell(x_train_projected)
        x_train_proj_mat = reshape(cell2mat(x_train_projected)', [ncomps^2, n_train_obs])';
        x_test_proj_mat = reshape(cell2mat(x_test_projected)', [ncomps^2, n_test_obs])';
    else
        x_train_proj_mat = reshape(x_train_projected, [ncomps^2, n_train_obs])';
        x_test_proj_mat = reshape(x_test_projected, [ncomps^2, n_test_obs])';
    end
else
    x_train_proj_mat = get_diag_elems(x_train_projected, n_train_obs, ncomps);
    x_test_proj_mat = get_diag_elems(x_test_projected, n_test_obs, ncomps);
    
end

mnr_mod = mnrfit(x_train_proj_mat, y_train);
predictions = mnrval(mnr_mod, x_test_proj_mat);
end

function diag_elems = get_diag_elems(x_projected, nobs, ncomps)
diag_elems = NaN([nobs, ncomps]);
for isample = 1:nobs
    if iscell(x_projected)
        diag_elems(isample, :) = diag(x_projected{isample});
    else
        diag_elems(isample, :) = diag(squeeze(x_projected(:, :, isample)));
    end
end

end