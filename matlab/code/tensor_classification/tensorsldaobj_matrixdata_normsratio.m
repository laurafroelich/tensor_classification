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

%trUtBU = sum(squeeze(sum(sum(between_classes_projected.^2, 1), 2)).*nis');
trUtBU = sum(sum(reshape(between_classes_projected.^2, [prod(Ks), length(nis)])).*nis);

% trace(U' * within_class_scatter * U)
% trUtWU = sum(sum(sum(between_obs_projected.^2)));
trUtWU = between_obs_projected.^2;
trUtWU = sum(trUtWU(:));


for imode = 1:nmodes
    temp_projected_classes = tmult(classmeandiffstensor,Us{2}',2);
    temp_projected_obs = tmult(observationdiffstensor,Us{2}',2);
    
    Ap2 = tmult(classmeandiffstensor,Us{1}',1);
    Bp2 = tmult(observationdiffstensor,Us{1}',1);
    
    S1=zeros([mode_sizes(1), Ks(1)]);
    S2=zeros([mode_sizes(2), Ks(2)]);
    for c=1:nclasses
        S1=S1+temp_projected_classes(:,:,c)*between_classes_projected(:,:,c)'*nis(c);
        S2=S2+Ap2(:,:,c)'*between_classes_projected(:,:,c)*nis(c);
    end
    
    T1=zeros([mode_sizes(1), Ks(1)]);
    T2=zeros([mode_sizes(2), Ks(2)]);
    for o=1:nobs
        T1=T1+temp_projected_obs(:,:,o)*between_obs_projected(:,:,o)';
        T2=T2+Bp2(:,:,o)'*between_obs_projected(:,:,o);
    end
    
    
    G1 = -(trUtWU*2*S1-trUtBU*2*T1)/trUtWU^2;
    G2 = -(trUtWU*2*S2-trUtBU*2*T2)/trUtWU^2;
    G.U1 = G1;
    G.U2 = G2;
end

F = -trUtBU/trUtWU;

end