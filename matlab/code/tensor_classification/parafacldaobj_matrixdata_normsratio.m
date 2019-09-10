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
trUtAU = sum(sum(reshape(between_classes_projected(all_inds).^2, [K, length(nis)])).*nis);
% for two dimensions:
% sum(squeeze(sum(sum(mAAp1.*between_classes_projected.^2, 1), 2)).*nis');

all_inds = get_level_wise_diag_inds(K, nobs);
trUtBU = sum(sum(reshape(between_obs_projected(all_inds).^2, [K, nobs])));

S_inds.type = '()';
subs = repmat({':'}, 1, nmodes);

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
    
    repmat_vector = ones(1, nmodes);
    repmat_vector(1) = mode_sizes(imode);
    
    diagonal_indices = get_level_wise_diag_inds(K, nmodes-1);
    
    T = zeros(size(Us{imode}));
    for o=1:nobs
        S_inds.subs = [subs, o];
        temp_slice = subsref(between_obs_projected, S_inds);
        diagonal_elements = temp_slice(diagonal_indices);
        T = T + ...
            subsref(temp_projected_obs, S_inds).*...
            repmat(diagonal_elements', repmat_vector);
        %diag(diagonal_elements);
        %repmat(reshape(diagonal_elements, [1, 3, 1]), [mode_sizes(imode), 1, 3]);
            %repmat(reshape(diagonal_elements, [1, 1, 3]), [mode_sizes(imode), 3, 1]); % * ... %temp_projected_obs(:,:,o) * ...
            %diag(diag(...
            %subsref(between_obs_projected, S_inds) ... %between_obs_projected(:,:,o)...
            %));
    end
    
    S = zeros(size(Us{imode}));
    for c=1:nclasses
        S_inds.subs = [subs, c];
        temp_slice = subsref(between_classes_projected, S_inds);
        diagonal_elements = temp_slice(diagonal_indices);
        
        S = S + ...
            subsref(temp_projected_classes, S_inds).*...
            repmat(diagonal_elements', repmat_vector)...         %diag(diagonal_elements)...
         *nis(c);   
        %repmat(reshape(diagonal_elements, [1, 3, 1]), [mode_sizes(imode), 1, 3]);
            %repmat(reshape(diagonal_elements, [1, 1, 3]), [mode_sizes(imode), 3, 1]);
        %S = S + ...
        %    subsref(temp_projected_classes, S_inds) * ... %temp_projected_classes(:,:,c) * ...
        %    diag(diag(...
        %    subsref(between_classes_projected, S_inds) ... %between_classes_projected(:,:,c)...
        %    ))*nis(c);
    end
    
    G.(['U', num2str(imode)]) = -(trUtBU*2*S-trUtAU*2*T)/trUtBU^2;
end

F = -trUtAU/trUtBU;

end