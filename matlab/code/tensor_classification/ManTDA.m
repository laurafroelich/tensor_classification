function [Us, outputs, Ys] = ManTDA(Xs, classes, varargin)
% [Us, outsmode1, outsmode2, outerwhile, Ys] = ManTDA(Xs, classes, varargin)
% classes: vector containing class labels. Classes must be sequential
% numbers starting from one.
% varargin{1}: lowerdims
% varargin{2}: Us
% varargin{3}: opts
% varargin{4}: stoppingcrit
% varargin{5}: maxits
% varargin{6}: opts
% varargin{7}: manopt method
%

%rng(0)
if isa(Xs, 'cell')
    Xs = cell_array_to_nd_array(Xs);
end

[nsamples, sizeX, nmodes] = get_sizes(Xs, 1); % observations assumed to run along first mode

Fdifftol = 1e-10;
Udifftol = 1e-12;
maxits = 1000;
Us = [];
optmeth = 'ManOpt';

lowerdims = sizeX;


usestoppingcrit = true;
opts = [];

if ~isempty(varargin)
    if ~isempty(varargin{1})
        lowerdims = varargin{1};
    end
    
    if length(varargin)>=2 && ~isempty(varargin{2})
        Us = varargin{2};
    end
    
    if length(varargin)>=3 && ~isempty(varargin{3})
        opts = varargin{3};
    end
    
    if length(varargin)>=4 && ~isempty(varargin{4})
        usestoppingcrit = varargin{4};
    end
    
    if length(varargin)>=5 && ~isempty(varargin{5})
        maxits = varargin{5};
    end
    
    if length(varargin)>=6 && ~isempty(varargin{6})
        store = varargin{6};
    end
    
    if length(varargin)>=7 && ~isempty(varargin{7})
        if ismember(varargin{7}, {'ManOpt', 'bo13', 'wen12'})
        optmeth = varargin{7};
        else
           warning('ManTDA.m: chosen manifold optimiser is not an option. Using ManOpt.') 
        end
    end
end


if isempty(Us)
    Us = cell(1, nmodes);
    for kmode = 1:nmodes
        Us{kmode} = orth(randn(sizeX(kmode), lowerdims(kmode)));
    end
end


Ks = [];
for imode = 1:nmodes
   U.(['U', num2str(imode)]) = Us{imode};
   Ks(imode) = size(Us{imode}, 2);
end

U.U1 = Us{1};
U.U2 = Us{2};
[N, K1] = size(Us{1});
[M, K2] = size(Us{2});
%U = zeros(N+M, K1+K2);
%U(1:N, 1:K1) = Us{1};
%U((N+1):end, (K1+1):end) = Us{2};

[cmean_m_xmeans, xi_m_cmeans, nis] = classbased_differences(Xs, classes);


% calculate Rw and Rb
[~, ~, Rw, Rb] = tensorsldaobj_matrixdata(U,...
    cmean_m_xmeans, xi_m_cmeans, nis, Ks);

opts.intialtau = -1;
opts.mxitr = maxits;
opts.record = 1;

options.maxiter = maxits;


switch optmeth
    case 'ManOpt'
        for imode = 1:nmodes
            tuple.(['U', num2str(imode)]) = stiefelfactory(size(Us{imode}, 1), size(Us{imode}, 2));
        end
        manifold = productmanifold(tuple);
        problem.M = manifold;
        
        % Define the problem cost function and its Euclidean gradient.
        problem.cost  = @(U, store) mycost(U, store,...
            cmean_m_xmeans, xi_m_cmeans, nis,...
            Ks, Rw, Rb);
        
        problem.egrad = @(U, store) mygrad(U, store,...
            cmean_m_xmeans, xi_m_cmeans, nis,...
            Ks, Rw, Rb);
        
        
        % Solve.
        [U, ~, outs] = conjugategradient(problem, U, options);
        fvals = cell2mat({outs.cost});
        
    case 'bo13'
        [U, outs] = OptStiManAFBB_myvariant(U,...
            @tensorsldaobj_matrixdata, opts, cmean_m_xmeans, xi_m_cmeans,...
            nis, Ks, Rw, Rb);
        fvals = outs.FArray;
        
    case 'wen12'
        [U, outs] = OptStiefelGBB_myvariant(U,...
            @tensorsldaobj_matrixdata, opts, cmean_m_xmeans, xi_m_cmeans,...
            nis, Ks, Rw, Rb);
        fvals = outs.fvals;
        
    otherwise
        error(['ManTDA.m: unknown optimisation method (' optmeth ').'])
end

outputs.fvals = fvals;
outputs.outs = outs;

for imode = 1:nmodes
   Us{imode} = U.(['U', num2str(imode)]);
end

if nargout >= 5
    Ys = cell(1, nsamples);
    for isample = 1:nsamples
        Yi = Xs{isample};
        for jmode = 1:nmodes % do multiplication between Xi tensor and U matrices
            origdims = size(Yi);
            Yi = Us{jmode}'*matricizing(Yi, jmode);
            newdims  = origdims;
            newdims(jmode) = size(Us{jmode}, 2);
            Yi = unmatricizing(Yi, jmode, newdims);
        end
        Ys{isample} = Yi;
    end
end
end

function [F, store] = mycost(x, store, classmeandiffs, observationdiffs,...
    nis, Ks, Rw, Rb)
[F, ~, ~, ~, store]...
    = tensorsldaobj_matrixdata(x,...
    classmeandiffs, observationdiffs, nis, Ks, ...
    Rw, Rb, store);
end


function [G, store] = mygrad(x, store, classmeandiffs, observationdiffs,...
    nis, Ks, Rw, Rb)
[~, G, ~, ~, store]...
    = tensorsldaobj_matrixdata(x,...
    classmeandiffs, observationdiffs, nis, Ks, ...
    Rw, Rb, store);
end

