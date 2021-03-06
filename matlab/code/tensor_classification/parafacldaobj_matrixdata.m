function [F, G, Rw, Rb, store] = parafacldaobj_matrixdata(U,...
    classmeandiffs, observationdiffs, nis, Ks, ...
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
% K:               Number of components in each mode.
%
% Rw and Rb:        Optional input. Pre-calculated matrices.

K = Ks(1);

if ~all(Ks == K)
   warning(['parafacldaobj_matrixdata.m: for PARAFAC, all projection ', ...
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

if ~storeexists || ~isfield(store, 'Rw') || ~isfield(store, 'Rb')
    
    if nargin < 7
        [~, ~, nmodes] = get_sizes(classmeandiffs, 1); % observations assumed to run along first mode
        
        permute_vector_move_obs_to_last_mode = [2:(nmodes+1), 1]; % move observations to run along last mode
        classmeandiffstensor = permute(classmeandiffs, ...
            permute_vector_move_obs_to_last_mode);
        observationdiffstensor = permute(observationdiffs, ...
            permute_vector_move_obs_to_last_mode);
        
        Rw = observationdiffstensor;
        Rb = classwise_scalar_multiply(classmeandiffstensor, sqrt(nis));
    end
    store.Rw = Rw;
    store.Rb = Rb;
end
Rw = store.Rw;
Rb = store.Rb;

Us = cell(1, nmodes);
for imode = 1:nmodes
    Us{imode} = U.(['U', num2str(imode)]);
end

permute_vector = [nmodes:-1:1 nmodes+1];
Rw_reversed_mode_order = permute(Rw, permute_vector);
Rb_reversed_mode_order = permute(Rb, permute_vector);

if ~isfield(store, 'Q')
    % calculate all scalars in the Kronecker product of all projection
    % matrices
    Q = reshape(Us{nmodes}, [mode_sizes(nmodes) 1 K]);
    for imode = (nmodes-1):-1:1
        A = reshape(Us{imode}, [1 mode_sizes(imode) K]);
        Q = reshape(bsxfun(@times,A,Q), [prod(mode_sizes(end:-1:imode)) 1 K]);
    end
    Q = reshape(Q, [size(Q, 1) K]);
    store.Q = Q;
else
    Q = store.Q;
end

if ~isfield(store, 'QtRw')
    QtRw = Q' * reshape(Rw_reversed_mode_order, mode_size_product, nobs);
    store.QtRw = QtRw;
else
    QtRw = store.QtRw;
end

if ~isfield(store, 'QtRb')
    QtRb = Q' * reshape(Rb_reversed_mode_order, mode_size_product, nclasses);
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

if ~isfield(store, 'QtWQinvQtBQ')
    QtWQinvQtBQ = (QtWQ)\(QtBQ);
    store.QtWQinvQtBQ = QtWQinvQtBQ;
else
    QtWQinvQtBQ = store.QtWQinvQtBQ;
end

F = trace(QtWQinvQtBQ);

if ~isfield(store, 'FdwrtQ')
    if mode_size_product < nobs % perform multiplication in fastest order
        FdwrtQ = (-2*(reshape(Rw_reversed_mode_order, ...
            mode_size_product, nobs)*QtRw')*QtWQinvQtBQ + ...
            2*reshape(Rb_reversed_mode_order, ...
            mode_size_product, nclasses)*QtRb')/QtWQ;
    else
        FdwrtQ = (-2*reshape(Rw_reversed_mode_order, ...
            mode_size_product, nobs)*(QtRw'*QtWQinvQtBQ) + ...
            2*reshape(Rb_reversed_mode_order, ...
            mode_size_product, nclasses)*QtRb')/QtWQ;
    end
    store.FdwrtQ = FdwrtQ;
else
    FdwrtQ = store.FdwrtQ;
end

for imode=1:nmodes
    permute_vector = [imode nmodes:-1:(imode+1) (imode-1):-1:1 nmodes+1];
    
    cost_derivative_new_size = [mode_sizes(end:-1:1), K];
    cost_derivative_2d_shape = mode_sizes;
    cost_derivative_2d_shape(imode) = cost_derivative_2d_shape(imode)*K;
    
    TTT=reshape(...
        permute(...
        reshape(FdwrtQ, cost_derivative_new_size),...
        [nmodes:-1:1 nmodes+1]),...
        cost_derivative_2d_shape);
    
    TTT3d = reshape(TTT, [mode_sizes, K]);
    G_imode = zeros(mode_sizes(imode), K);
    modes_to_multiply = setdiff(1:nmodes, imode);
    temp_tmult = TTT3d;
    for mode = modes_to_multiply
        temp_tmult = tmult(temp_tmult, Us{mode}', mode);
    end
    temp_tmult = permute(temp_tmult, permute_vector);
    
    S.type = '()';
    for ilowerdim = 1:K
        S.subs = [{':'} repmat([{ilowerdim}], 1, nmodes)];
        G_imode(:, ilowerdim) = subsref(temp_tmult,S);
    end
    G.(['U', num2str(imode)]) = -G_imode;
end
F = -F;
end