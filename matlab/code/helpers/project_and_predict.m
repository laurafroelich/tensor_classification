function predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us, varargin)

if isempty(varargin)
    parafac_structure = false;
else
    parafac_structure = varargin{1};
end

nmodes = length(Us);
x_train_projected = tensor_projection(x_train, Us);
x_test_projected = tensor_projection(x_test, Us);

n_train_obs = length(y_train);
n_test_obs = length(y_test);

if ~parafac_structure
    if iscell(x_train_projected)
        permute_vector = [nmodes+1, 2:nmodes, 1];
        
        x_train_proj_mat = cell_array_to_nd_array(x_train_projected); %reshape(cell2mat(x_train_projected)', [ncomps^2, n_train_obs])';
        x_test_proj_mat = cell_array_to_nd_array(x_test_projected); %reshape(cell2mat(x_test_projected)', [ncomps^2, n_test_obs])';
        x_train_proj_mat = reshape(permute(x_train_proj_mat, permute_vector),...
            [prod(ncomps), n_train_obs])';
        x_test_proj_mat = reshape(permute(x_test_proj_mat, permute_vector),...
            [prod(ncomps), n_test_obs])';
        %x_train_proj_mat = reshape(cell2mat(x_train_projected)', [ncomps^2, n_train_obs])';
        %x_test_proj_mat = reshape(cell2mat(x_test_projected)', [ncomps^2, n_test_obs])';
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
k=ncomps(1);
diag_elems = NaN([nobs, k]);

for isample = 1:nobs
    if iscell(x_projected)
        % use linear indices to get diagonal elements, as given here:
        % https://stackoverflow.com/questions/5598900/how-can-i-index-the-diagonals-of-a-3-d-matrix-in-matlab
        if k == length(ncomps)
        inds = get_diagonal_indices(k, length(ncomps));
        else
        inds = diag(cumsum([1:(k+1):k^2; k^2.*ones(k-1,k)]));
        end
        temp = x_projected{isample};
        diag_elems(isample, :) = temp(inds);
        %diag_elems(isample, :) = diag(x_projected{isample});
    else
        diag_elems(isample, :) = diag(squeeze(x_projected(:, :, isample)));
    end
end

end