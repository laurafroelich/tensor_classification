function [F, G, classmeandiffstensor, observationdiffstensor, store] = ...
    tensorsldaobj_matrixdata_normsratio(U,...
    classmeandiffs, observationdiffs, nis, Ks,...
    classmeandiffstensor, observationdiffstensor, store)
% [F, G, Ap, Bp] = tensorsldaobj_matrixdata_normsratio(U,...
%    classmeandiffs, observationdiffs, nis, K1, K2,...
%    Ap, Bp)
%
% U:
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
% K1:               Number of factors in mode 1.
%
% K2:               Number of factors in mode 2.
% classmeandiffstensor and observationdiffstensor:        Optional input. Pre-calculated matrices.

[nclasses, mode_sizes, nmodes] = get_sizes(classmeandiffs, 1); % observations assumed to run along first mode
nobs = sum(nis);

Us = cell(1, nmodes);
for imode = 1:nmodes
    Us{imode} = U.(['U', num2str(imode)]);
end

if nargin <=6
    permute_vector = [2:(nmodes+1), 1];
    classmeandiffstensor = permute(classmeandiffs, permute_vector);
    observationdiffstensor = permute(observationdiffs, permute_vector);
end


between_classes_projected = classmeandiffstensor;
between_obs_projected = observationdiffstensor;
for imode = 1:nmodes
    between_classes_projected = tmult(between_classes_projected, Us{imode}', imode);
    between_obs_projected = tmult(between_obs_projected, Us{imode}', imode);
end

% variable name: trace(U' * between_class_scatter * U)
trUtBU = sum(sum(reshape(between_classes_projected.^2, [prod(Ks), length(nis)])).*nis);

% variable name: trace(U' * within_class_scatter * U)
trUtWU = between_obs_projected.^2;
trUtWU = sum(trUtWU(:));


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
    between_obs_projected_permuted = permute(between_obs_projected, permute_vector2);
    
    temp_projected_classes = permute(temp_projected_classes, permute_vector2);
    between_classes_projected_permuted = permute(between_classes_projected, permute_vector2);
    between_classes_projected_permuted = between_classes_projected_permuted .* ...
        reshape(repelem(nis, prod(Ks)), [Ks, nclasses]);
    
    S = reshape(temp_projected_classes, [mode_sizes(imode), prod([nclasses, Ks(othermodes)])]) * ...
        reshape(between_classes_projected_permuted, [Ks(imode), prod([nclasses, Ks(othermodes)])])';
    
    T = reshape(temp_projected_obs, [mode_sizes(imode), prod([nobs, Ks(othermodes)])]) * ...
        reshape(between_obs_projected_permuted, [Ks(imode), prod([nobs, Ks(othermodes)])])';
    
    G.(['U', num2str(imode)]) = -(trUtWU*2*S-trUtBU*2*T)/trUtWU^2;
end

F = -trUtBU/trUtWU;

end