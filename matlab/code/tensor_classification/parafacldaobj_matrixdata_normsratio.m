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

if nargin < 9
    storeexists = false;
else
    storeexists = true;
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

if ~storeexists || ~isfield(store, 'Ap1')
    Ap1 = tmult(classmeandiffstensor,Us{2}',2);
    store.Ap1 = Ap1;
else
    Ap1 = store.Ap1;
end

if ~isfield(store, 'Bp1')
    Bp1 = tmult(observationdiffstensor,Us{2}',2);
    store.Bp1 = Bp1;
else
    Bp1 = store.Bp1;
end

if ~isfield(store, 'Ap2')
    Ap2 = tmult(classmeandiffstensor,Us{1}',1);
    store.Ap2 = Ap2;
else
    Ap2 = store.Ap2;
end

if ~isfield(store, 'Bp2')
    Bp2 = tmult(observationdiffstensor,Us{1}',1);
    store.Bp2 = Bp2;
else
    Bp2 = store.Bp2;
end

if ~isfield(store, 'AAp1')
    AAp1 = tmult(Ap1,Us{1}',1);
    store.AAp1 = AAp1;
else
    AAp1 = store.AAp1;
end

if ~isfield(store, 'BBp1')
    BBp1 = tmult(Bp1,Us{1}',1);
    store.BBp1 = BBp1;
else
    BBp1 = store.BBp1;
end

if ~isfield(store, 'AAp2')
    AAp2 = tmult(Ap2,Us{2}',2);
    store.AAp2 = AAp2;
else
    AAp2 = store.AAp2;
end

if ~isfield(store, 'BBp2')
    BBp2 = tmult(Bp2,Us{2}',2);
    store.BBp2 = BBp2;
else
    BBp2 = store.BBp2;
end

if ~isfield(store, 'mAAp1')
    mAAp1=repmat(eye(size(AAp1,1)),[1 1 size(AAp1,3)]);
    store.mAAp1 = mAAp1;
else
    mAAp1 = store.mAAp1;
end

if ~isfield(store, 'mBBp1')
    mBBp1=repmat(eye(size(BBp1,1)),[1 1 size(BBp1,3)]);
    store.mBBp1 = mBBp1;
else
    mBBp1 = store.mBBp1;
end

if ~isfield(store, 'trUtAU')
    trUtAU = sum(squeeze(sum(sum(mAAp1.*AAp1.^2, 1), 2)).*nis');
    store.trUtAU = trUtAU;
else
    trUtAU = store.trUtAU;
end

if ~isfield(store, 'trUtBU')
    trUtBU = sum(sum(sum(mBBp1.*BBp1.^2)));
    store.trUtBU = trUtBU;
else
    trUtBU = store.trUtBU;
end


if ~isfield(store, 'Ss')
    Ss = cell(1, nmodes);
    for imode = 1:nmodes
        Ss{imode} = zeros(size(Us{imode}));
    end
    
    for c=1:size(Ap1,3)
        Ss{1}=Ss{1}+Ap1(:,:,c)*diag(diag(AAp1(:,:,c)))*nis(c);
        Ss{2}=Ss{2}+Ap2(:,:,c)'*diag(diag(AAp2(:,:,c)))*nis(c);
    end
    store.Ss = Ss;
else
    Ss = store.Ss;
end

if ~isfield(store, 'Ts')
    
    Ts = cell(1, nmodes);
    for imode = 1:nmodes
        Ts{imode} = zeros(size(Us{imode}));
    end
    for o=1:size(Bp1,3)
        Ts{1}=Ts{1}+Bp1(:,:,o)*diag(diag(BBp1(:,:,o)))';
        Ts{2}=Ts{2}+Bp2(:,:,o)'*diag(diag(BBp2(:,:,o)));
    end
    store.Ts = Ts;
else
    Ts = store.Ts;
end

for imode = 1:nmodes
    G.(['U', num2str(imode)]) = -(trUtBU*2*Ss{imode}-trUtAU*2*Ts{imode})/trUtBU^2;
end

F = -trUtAU/trUtBU;

end