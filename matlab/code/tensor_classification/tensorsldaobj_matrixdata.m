function [F, G, Rw, Rb, store]...
    = tensorsldaobj_matrixdata(U,...
    classmeandiffs, observationdiffs, nis, Ks, ...
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

[nclasses, mode_sizes, nmodes] = get_sizes(classmeandiffs, 1); % observations assumed to run along first mode
mode_size_product = prod(mode_sizes);

if nargin < 9
    storeexists = false;
else
    storeexists = true;
end

if ~storeexists || ~isfield(store, 'Rw') || ~isfield(store, 'Rb')
    
    if nargin < 7
        
        permute_vector_move_obs_to_last_mode = [2:(nmodes+1), 1]; % move observations to run along last mode
        classmeandiffstensor = permute(classmeandiffs, permute_vector_move_obs_to_last_mode);
        observationdiffstensor = permute(observationdiffs, permute_vector_move_obs_to_last_mode);
        
        Rw =observationdiffstensor;
        
        Rb = classwise_scalar_multiply(classmeandiffstensor, sqrt(nis));
    end
    store.Rw = Rw;
    store.Rb = Rb;
end

Rw = store.Rw;
Rb = store.Rb;

Rwsize = size(Rw);
nobs = Rwsize(end);

Us = cell(1, nmodes);
for imode = 1:nmodes
    Us{imode} = U.(['U', num2str(imode)]);
end

permute_vector = [nmodes:-1:1 nmodes+1];

QtRw = Rw;
for imode = 1:nmodes
    QtRw = tmult(QtRw, Us{imode}', imode);
end

QtRw = permute(QtRw, permute_vector);
QtRw = reshape(QtRw, prod(Ks), nobs);

QtRb = Rb;
for imode = 1:nmodes
    QtRb = tmult(QtRb, Us{imode}', imode);
end
QtRb = permute(QtRb, permute_vector);
QtRb = reshape(QtRb, prod(Ks), nclasses);

QtWQ = QtRw*QtRw';

QtBQ = QtRb*QtRb';

QtWQinvQtBQ = (QtWQ)\(QtBQ);

F = trace(QtWQinvQtBQ);

if mode_size_product < nobs % perform multiplication in fastest order
    FdwrtQ = (-2*(reshape(permute(Rw, permute_vector), ...
        mode_size_product, nobs)*QtRw')*QtWQinvQtBQ + ...
        2*reshape(permute(Rb, permute_vector), ...
        mode_size_product, nclasses)*QtRb')/QtWQ;
else
    FdwrtQ = (-2*reshape(permute(Rw, permute_vector), ...
        mode_size_product, nobs)*(QtRw'*QtWQinvQtBQ) + ...
        2*reshape(permute(Rb, permute_vector), ...
        mode_size_product, nclasses)*QtRb')/QtWQ;
end

% Use this permutation vector to change the order of modes such that
% the mode that was highest is moved to lowest and vice versa. Also,
% follow each mode dimensions by its corresponding representation in
% the lower projected dimension.
permute_vector_mode_followed_by_lower_dims = [];
for temp_mode = 1:nmodes
    permute_vector_mode_followed_by_lower_dims = ...
        [temp_mode, temp_mode+nmodes, ...
        permute_vector_mode_followed_by_lower_dims];
end

for imode=1:nmodes
    
    cost_derivative_new_size = [mode_sizes(end:-1:1), Ks(end:-1:1)];
    cost_derivative_2d_shape = mode_sizes.*Ks;
    
    TTT=reshape(...
        permute(...
        reshape(FdwrtQ, cost_derivative_new_size),...
        permute_vector_mode_followed_by_lower_dims),...
        cost_derivative_2d_shape);
    
    modes_to_multiply = setdiff(1:nmodes, imode);
    temp_tmult = TTT;
    for mode = modes_to_multiply
        U = Us{mode};
        temp_tmult = tmult(temp_tmult, U(:)', mode);
    end
    G.(['U', num2str(imode)]) = -reshape(temp_tmult, mode_sizes(imode), Ks(imode));
end
F = -F;

end


