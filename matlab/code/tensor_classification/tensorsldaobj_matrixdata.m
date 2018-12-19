function [F, G, Rw, Rb, store]...
    = tensorsldaobj_matrixdata(U,...
    classmeandiffs, observationdiffs, nis, K1, K2, ...
    Rw, Rb, store)
%[F, G, Rw, Rb, store]...
%    = tensorsldaobj_matrixdata(U,...
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

if ~isfield(store, 'QtRw')
    QtRw_mm=tmult(tmult(Rw,U1',1),U2',2);
    QtRw=reshape(permute(QtRw_mm,[2 1 3]),[K2*K1,nobs]);
    %mynumb=25; cond(reshape(permute(QtRw_mm(:,:,1:mynumb),[2 1 3]),[L*K,mynumb]))
    store.QtRw = QtRw;
else
    QtRw = store.QtRw;
end

if ~isfield(store, 'QtRb')
    QtRb_mm=tmult(tmult(Rb,U1',1),U2',2);
    QtRb=reshape(permute(QtRb_mm,[2 1 3]),[K2*K1,nclasses]);
    %mynumb=1; cond(reshape(permute(QtRb_mm(:,:,1:mynumb),[2 1 3]),[L*K,mynumb]))
    
    store.QtRb = QtRb;
else
    QtRb = store.QtRb;
end

if ~isfield(store, 'QtWQ')
    QtWQ = QtRw*QtRw';
    store.QtWQ = QtWQ;
else
    QtWQ = store.QtWQ;
end

if ~isfield(store, 'QtBQ')
    QtBQ = QtRb*QtRb';
    store.QtBQ = QtBQ;
else
    QtBQ = store.QtBQ;
end

if ~isfield(store, 'QtWQinvQtBQ')
    QtWQinvQtBQ = (QtWQ)\(QtBQ);
    store.QtWQinvQtBQ = QtWQinvQtBQ;
else
    QtWQinvQtBQ = store.QtWQinvQtBQ;
end

F = trace(QtWQinvQtBQ);

% inefficient calculations
if false
    Q = kron(U1, U2);
    Qt= Q';
    W = zeros(I*J, I*J);
    for iobs = 1:nobs
        curobs = observationdiffstensor(:,:,iobs);
        tempw = reshape(curobs', I*J, 1);
        W = W + tempw*tempw';
    end
    B = zeros(I*J, I*J);
    for iclass = 1:nclasses
        curclass = classmeandiffstensor(:,:,iclass);
        tempb = reshape(curclass', I*J, 1);
        B = B + nis(iclass)*tempb*tempb';
    end
    QtWQx = Qt*W*Q;
    QtBQx = Qt*B*Q;
    QtWQinvQtBQx = (QtWQx)\(QtBQx);
    Fx = trace(QtWQinvQtBQx);
end

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

if ~isfield(store, 'TTT')
    TTT=reshape(permute(reshape(FdwrtQ, M, N, K2, K1), [2 4 1 3]),[N*K1 M*K2]);
    store.TTT = TTT;
else
    TTT = store.TTT;
end

G1 =reshape(TTT*U2(:),[N K1]);
G2=reshape(TTT'*U1(:),[M K2]);

F = -F;
%G = zeros(N+M, K1+K2);
%G(1:N, 1:K1) = G1;
%G((N+1):end, (K1+1):end) = G2;
%G=-G;
G.U1 = -G1;
G.U2 = -G2;
end


