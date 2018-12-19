function [F, G, Rw, Rb, store] = parafacldaobj_matrixdata(U,...
    classmeandiffs, observationdiffs, nis, K1, K2, ...
    Rw, Rb, store)
%[F, G, Rw, Rb, store]...
%    = parafacldaobj_matrixdata(U,...
%    classmeandiffs, observationdiffs, nis, K1, K2, ...
%    Rw, Rb, store)

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
        obsexample = classmeandiffs{1};
        sizeobs = size(obsexample);
        I = sizeobs(1);
        J = sizeobs(2);
        nclasses = length(classmeandiffs);
        nobs = length(observationdiffs);
        classmeandiffstensor = reshape(cell2mat(classmeandiffs), ...
            I, J, nclasses);
        observationdiffstensor = reshape(cell2mat(observationdiffs), ...
            I, J, nobs);
        
        Rw =observationdiffstensor;
        Rb = classmeandiffstensor.*permute(repmat(sqrt(nis), I,1,J), [1 3 2]);
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
%U1=U(1:N, 1:K1);
%U2=U((N+1):end, (K1+1):end);
U1 = U.U1;
U2 = U.U2;

%QtRw_mm=tmult(tmult(Rw,U1',1),U2',2);
%QtRw=reshape(permute(QtRw_mm,[2 1 3]),[L*K,nobs]);

%QtRb_mm=tmult(tmult(Rb,U1',1),U2',2);
%QtRb=reshape(permute(QtRb_mm,[2 1 3]),[L*K,nclasses]);

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
    %QtWQ = QtRw*QtRw';
    QtWQ = diag(diag(QtRw*QtRw'));
    store.QtWQ = QtWQ;
else
    QtWQ = store.QtWQ;
end
if ~isfield(store, 'QtBQ')
    %QtBQ = QtRb*QtRb';
    QtBQ = diag(diag(QtRb*QtRb'));
    store.QtBQ = QtBQ;
else
    QtBQ = store.QtBQ;
end

if false % alternative, less efficient way to calculate QtWQ
    QtRw_mmalt=zeros(K, nobs);%tmult(tmult(Rw,U1(:, 1)',1),U2(:,1)',2);
    for curcomp = 1:K
       QtRw_mmalt(curcomp, :) = squeeze(tmult(tmult(Rw,U1(:, curcomp)',1),U2(:,curcomp)',2))';
    end
    %QtRwalt=reshape(permute(QtRw_mmalt,[2 1 3]),[K^2,nobs]);
    sum(QtRw_mmalt.^2,2)
    %QtWQalt = QtRwalt*QtRwalt';
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
%Gtemp = reshape(TTT*U2, [N, L, K]);
%arrayfun(@(x)(Gtemp(:, x, x)), 1:K, 'UniformOutput', false)
G1 = cell2mat(arrayfun(@(x)(TTT3d(:,:,x)*U2(:,x)), 1:K, 'UniformOutput', false));

%G =reshape(TTT*U2,[N K]); % we need the first N numbers from column 1, the next N numbers from column 2, etc.


TTT=reshape(permute(reshape(FdwrtQ, M, N, K, 1), [2 4 1 3]),[N M*K]);
TTT3d = permute(reshape(TTT, [N, M, K]), [2 1 3]);
G2 = cell2mat(arrayfun(@(x)(TTT3d(:,:,x)*U1(:,x)), 1:K, 'UniformOutput', false));



F = -F;
%G = zeros(N+M, K1+K2);
%G(1:N, 1:K1) = G1;
%G((N+1):end, (K1+1):end) = G2;
%G=-G;
G.U1 = -G1;
G.U2 = -G2;
end











