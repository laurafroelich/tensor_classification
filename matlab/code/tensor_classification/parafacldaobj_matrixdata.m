function [F, G, Rw, Rb, store] = parafacldaobj_matrixdata(U,...
    classmeandiffs, observationdiffs, nis, K1, K2, ...
    Rw, Rb, store)
%[F, G, Rw, Rb, store]...
%    = parafacldaobj_matrixdata(U,...
%    classmeandiffs, observationdiffs, nis, K1, K2, ...
%    Rw, Rb, store)

% U:
%
% classmeandiffs:   Array of class differences from mean of all data.
%
% observationdiffs: Array of I*J matrices holding X_i-M_c(i) for each
%                   observation X_i, where M_c(i) is the mean of X_i's
%                   class, c(i).
%
% nis:              Vector containing number of observations from each
%                   class.
%
% K1:               Number of components in mode 1.
%
% K2:               Number of components in mode 2.
%
% Rw and Rb:        Optional input. Pre-calculated matrices.

if K1~=K2
    error('parafacldaobj_matrixdata.m: number of components must be the same for the two modes.')
else
    K=K1;
end

if nargin < 9
    storeexists = false;
else
    storeexists = true;
end

if ~storeexists || ~isfield(store, 'Rw') || ~isfield(store, 'Rb')
    
    if nargin < 7
        [~, ~, nmodes] = get_sizes(classmeandiffs, 1); % observations assumed to run along first mode
        
        permute_vector = [2:(nmodes+1), 1]; % move observations to run along last mode
        classmeandiffstensor = permute(classmeandiffs, permute_vector);
        observationdiffstensor = permute(observationdiffs, permute_vector);
        
        Rw = observationdiffstensor;
        Rb = classwise_scalar_multiply(classmeandiffstensor, sqrt(nis));
    end
    store.Rw = Rw;
    store.Rb = Rb;
end
Rw = store.Rw;
Rb = store.Rb;

Rwsize = size(Rw);

nobs = Rwsize(end);
datadims = size(Rb);
nclasses = datadims(length(datadims));

N=datadims(1);
M=datadims(2);
U1 = U.U1;
U2 = U.U2;

if ~isfield(store, 'Q')
    Q = reshape(U2,[M 1 K]);
    A = reshape(U1,[1 N K]);
    Q = reshape(bsxfun(@times,A,Q),[N*M 1 K]);
    Q = reshape(Q,[size(Q,1) K]);
    store.Q = Q;
else
    Q = store.Q;
end
if ~isfield(store, 'QtRw')
    QtRw = Q' * reshape(permute(Rw, [2 1 3]), M*N, nobs);
    store.QtRw = QtRw;
else
    QtRw = store.QtRw;
end

if ~isfield(store, 'QtRb')
    QtRb = Q' * reshape(permute(Rb, [2 1 3]), M*N, nclasses);
    store.QtRb = QtRb;
else
    QtRb = store.QtRb;
end


if ~isfield(store, 'QtWQ')
    QtWQ = diag(diag(QtRw*QtRw'));
    store.QtWQ = QtWQ;
else
    QtWQ = store.QtWQ;
end
if ~isfield(store, 'QtBQ')
    QtBQ = diag(diag(QtRb*QtRb'));
    store.QtBQ = QtBQ;
else
    QtBQ = store.QtBQ;
end

if false % alternative, less efficient way to calculate QtWQ
    QtRw_mmalt=zeros(K, nobs);
    for curcomp = 1:K
       QtRw_mmalt(curcomp, :) = squeeze(tmult(tmult(Rw,U1(:, curcomp)',1),U2(:,curcomp)',2))';
    end
    sum(QtRw_mmalt.^2,2)
end

if ~isfield(store, 'QtWQinvQtBQ')
    QtWQinvQtBQ = (QtWQ)\(QtBQ);
    store.QtWQinvQtBQ = QtWQinvQtBQ;
else
    QtWQinvQtBQ = store.QtWQinvQtBQ;
end

F = trace(QtWQinvQtBQ);

if ~isfield(store, 'FdwrtQ')
    if M*N < nobs % perform multiplication in fastest order
        FdwrtQ = (-2*(reshape(permute(Rw,[2 1 3]), M*N, nobs)*QtRw')*QtWQinvQtBQ + 2*reshape(permute(Rb,[2 1 3]), M*N, nclasses)*QtRb')/QtWQ;
    else
        FdwrtQ = (-2*reshape(permute(Rw,[2 1 3]), M*N, nobs)*(QtRw'*QtWQinvQtBQ) + 2*reshape(permute(Rb,[2 1 3]), M*N, nclasses)*QtRb')/QtWQ;
    end
    store.FdwrtQ = FdwrtQ;
else
    FdwrtQ = store.FdwrtQ;
end



TTT=reshape(permute(reshape(FdwrtQ, M, N, 1, K), [2 4 1 3]),[N*K M]);
TTT3d = permute(reshape(TTT', [M, N, K]), [2 1 3]);
G1 = cell2mat(arrayfun(@(x)(TTT3d(:,:,x)*U2(:,x)), 1:K, 'UniformOutput', false));

TTT=reshape(permute(reshape(FdwrtQ, M, N, K, 1), [2 4 1 3]),[N M*K]);
TTT3d = permute(reshape(TTT, [N, M, K]), [2 1 3]);
G2 = cell2mat(arrayfun(@(x)(TTT3d(:,:,x)*U1(:,x)), 1:K, 'UniformOutput', false));

F = -F;
G.U1 = -G1;
G.U2 = -G2;
end











