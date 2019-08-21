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
mode_size_product = prod(mode_sizes);
nobs = sum(nis);

obsexample = classmeandiffs(1,:,:);
sizeobs = size(obsexample);

Us = cell(1, nmodes);
for imode = 1:nmodes
    Us{imode} = U.(['U', num2str(imode)]);
end

if nargin <=6
    permute_vector = [2:(length(sizeobs)), 1];
    classmeandiffstensor = permute(classmeandiffs, permute_vector);
    observationdiffstensor = permute(observationdiffs, permute_vector);
end

between_classes_projected = classmeandiffstensor;
between_obs_projected = observationdiffstensor;
for imode = 1:nmodes
    between_classes_projected = tmult(between_classes_projected, Us{imode}', imode);
    between_obs_projected = tmult(between_obs_projected, Us{imode}', imode);
end

mAAp1=repmat(eye(K),[1 1 nclasses]);

mBBp1=repmat(eye(K),[1 1 nobs]);

trUtAU = sum(squeeze(sum(sum(mAAp1.*between_classes_projected.^2, 1), 2)).*nis');

trUtBU = sum(sum(sum(mBBp1.*between_obs_projected.^2)));


for imode = 1:nmodes
    T = zeros(size(Us{imode}));
    
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
    
    for o=1:nobs
        T = T + ...
            temp_projected_obs(:,:,o) * ...
            diag(diag(between_obs_projected(:,:,o)));
    end
    
    S = zeros(size(Us{imode}));
    for c=1:nclasses
        S = S + ...
            temp_projected_classes(:,:,c) * ...
            diag(diag(between_classes_projected(:,:,c)))*nis(c);
    end
    
    G.(['U', num2str(imode)]) = -(trUtBU*2*S-trUtAU*2*T)/trUtBU^2;
end

F = -trUtAU/trUtBU;

end