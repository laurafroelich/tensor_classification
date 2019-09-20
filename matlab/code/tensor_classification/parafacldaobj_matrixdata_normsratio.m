function [F, G, classmeandiffstensor, observationdiffstensor, store] = ...
    parafacldaobj_matrixdata_normsratio(U,...
    classmeandiffs, observationdiffs, nis, Ks,...
    classmeandiffstensor, observationdiffstensor, store)
% [F, G, classmeandiffstensor, observationdiffstensor] = ...
%  parafacldaobj_matrixdata_normsratio(U,...
%    classmeandiffs, observationdiffs, nis, Ks,...
%    classmeandiffstensor, observationdiffstensor)
%
% Utooptimise:      If modetooptimise is 1, I*K1 matrix where K1 is the
%                   number of components for mode 1. Otherwise J*K2 matrix,
%                   where K2 is the number of components for mode 2.
%
% classmeandiffs:   Cell array of I*J matrices holding M_c-M for each
%                   class, where M_c is the mean of class c and M is the
%                   overall mean.
%
% observationdiffs: Cell array of I*J matrices holding X_i-M_c(i) for each
%                   observation X_i, where M_c(i) is the mean of X_i's
%                   class, c(i).
%
% nis:              Vector containing number of observations from each
%                   class.
%
% modetooptimise:   1 to optimise U for the first mode and 2 to optimise U
%                   for the second mode.
%
% otherU:           If modetooptimise is 1, J*K2 matrix where K2 is the
%                   number of components for mode 2. Otherwise I*K1 matrix,
%                   where K1 is the number of components for mode 1.
%
% Ap and Bp:        Optional input. Pre-calculated matrices.

K = Ks(1);

if ~all(Ks == K)
    warning(['parafacldaobj_matrixdata_normsratio.m: for PARAFAC, all projection ', ...
        'dimensions must be the same size, but ', num2str(Ks), ...,
        ' were given as sizes of lower dimensions'])
end

[nclasses, mode_sizes, nmodes] = get_sizes(classmeandiffs, 1); % observations assumed to run along first mode
nobs = sum(nis);

Us = cell(1, nmodes);
for imode = 1:nmodes
    Us{imode} = U.(['U', num2str(imode)]);
end

if nargin <=6
    permute_vector = [2:(nmodes + 1), 1];
    classmeandiffstensor = permute(classmeandiffs, permute_vector);
    observationdiffstensor = permute(observationdiffs, permute_vector);
end

between_classes_projected = classmeandiffstensor;
between_obs_projected = observationdiffstensor;
for imode = 1:nmodes
    between_classes_projected = tmult(between_classes_projected, Us{imode}', imode);
    between_obs_projected = tmult(between_obs_projected, Us{imode}', imode);
end

all_inds = get_level_wise_diag_inds(K, length(nis));
trUtBU = sum(sum(reshape(between_classes_projected(all_inds).^2, [K, length(nis)])).*nis);

all_inds = get_level_wise_diag_inds(K, nobs);
trUtWU = sum(between_obs_projected(all_inds).^2);

for imode = 1:nmodes
    
    othermodes = setdiff(1:nmodes, imode);
    
    temp_projected_obs = observationdiffstensor;
    temp_projected_classes = classmeandiffstensor;
    for other_mode = othermodes
        temp_projected_obs = tmult(temp_projected_obs, Us{other_mode}', other_mode);
        temp_projected_classes = tmult(temp_projected_classes, Us{other_mode}', other_mode);
    end
    
    permute_vector2 = [1:nmodes, nmodes+1];
    permute_vector2(1) = imode;
    permute_vector2(imode) = 1;
    temp_projected_obs = permute(temp_projected_obs, permute_vector2);
    temp_projected_classes = permute(temp_projected_classes, permute_vector2);
    
    rep_vector = ones(1, nmodes+1);
    rep_vector(1) = mode_sizes(imode);
    rep_vector(2:(nmodes-1)) = K;
    rep_vector = mat2cell(rep_vector, 1, ones(1, numel(rep_vector)));
    
    diagonal_indices_obs = get_level_wise_diag_inds(K, nobs, nmodes);
    between_obs_projected_diagonal = reshape(...
        between_obs_projected(diagonal_indices_obs), [K, nobs]);
    
    scaled_projected_obs_nd = temp_projected_obs.*...
        reshape(repelem(between_obs_projected_diagonal, rep_vector{:}), ...
        [mode_sizes(imode), repmat(K, 1, nmodes-1), nobs]);
    scaled_projected_obs_nd_sum = sum(scaled_projected_obs_nd, nmodes+1);
    
    
    diagonal_indices_classes = get_level_wise_diag_inds(K, nclasses, nmodes);
    between_classes_projected_diagonal = reshape(...
        between_classes_projected(diagonal_indices_classes), [K, nclasses]);
    
    scaled_projected_classes_nd = temp_projected_classes.*...
        reshape(repelem(between_classes_projected_diagonal.*nis, rep_vector{:}), ...
        [mode_sizes(imode), repmat(K, 1, nmodes-1), nclasses]);
    scaled_projected_classes_nd_sum = sum(scaled_projected_classes_nd, nmodes+1);
    
    G.(['U', num2str(imode)]) = -(trUtWU*2*scaled_projected_classes_nd_sum...
        -trUtBU*2*scaled_projected_obs_nd_sum)/trUtWU^2;
    
end

F = -trUtBU/trUtWU;

end