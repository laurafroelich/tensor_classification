function [Us, iit, errs, objfuncvals, objfuncvals_traceratio, Ys] = CMDA(Xs, classes, varargin)
% [Us, iit, errs, objfuncvals, Ys] = CMDA(Xs, classes, varargin)
%
% Use CMDA to find projection matrices for each mode of the input data in Xs.
% Xs is assumed to be a cell array, with each entry containing an
% n-dimensional observation tensor.
%
% Mandatory input:
% Xs:           Cell array containing the observed tensors or array with 
%               observations along last dimension.
% classes:      Vector containing class labels. Classes must be sequential
%               numbers starting from one.
%
% Optional input:
% varargin{1}:  Scalar giving the maximal number of (outer) iterations.
%               Default: 100.
% varargin{2}:  Vector containing the number of components to estimate for
%               each mode. Default: size of observations.
% varargin{3}:  Boolean indicating whether or not convergence criteria
%               should be allowed to terminate iterations.
%               Default: true.
% varargin{4}:  Cell array of initial projection matrices or the string
%               'randinit' indicating initialisation with random
%               orthogonal matrices. Default: all-one matrices as
%               proposed in li14.
%
% Output:
% Us:           Cell array containing the projection matrices found.
% iit:          The number of outer iterations performed (an outer
%               iteration consists of one optimisation of each mode).
% errs:         The values from the stopping criterion proposed in li14,
%               one for each outer iteration.
% objfuncvals:  Values of the objective function that CMDA tries to
%               optimise. Values are given for all inner iterations. An
%               inner iteration consists of the optimisation of one mode.
% objfuncvals_traceratio:   Values of the objective function
%                           tr((U'WU)^(-1)U'BU). W is
%                           the within-class scatter matrix while B is the
%                           between-class scatter matrix.
% Ys:           Projections of the original data in Xs projected onto the
%               final projection matrices in Us.
%
% li14:
%   Q. Li and D. Schonfeld,
%   'Multilinear Discriminant Analysis for Higher-Order Tensor
%   Data Classification',
%   IEEE Transactions on Pattern Analysis and Machine Intelligence

%% read input and set parameters

if isa(Xs, 'cell')
    Xs = cell_array_to_nd_array(Xs);
end
    
[nobs, sizeX, nmodes] = get_sizes(Xs, 1); % observations assumed to run along first mode
tol=1e-6;

if length(varargin) >= 1 && ~isempty(varargin{1})
    Tmax = varargin{1};
else
    Tmax = 100;
end

if length(varargin)>=2 && ~isempty(varargin{2})
    lowerdims = varargin{2};
else
    lowerdims = sizeX;
end

if length(varargin)>=3 && ~isempty(varargin{3})
    usestoppingcrit = varargin{3};
else
    usestoppingcrit = true;
end

if length(varargin)>=4 && ~isempty(varargin{4})
    Us = varargin{4};
    if ischar(Us)
        switch Us
            case 'randinit'
                Us = cell(1, nmodes);
                for kmode = 1:nmodes
                    Us{kmode} = orth(randn(sizeX(kmode), lowerdims(kmode)));
                end
            otherwise
                warning(['CMDA.m: initialisation method not recognised, '...
                    'initialising with all-one matrices as proposed ' ...
                    'in li14 (see help for citation)'])
                % initialisation as proposed in li14
                Us = cell(1, nmodes);
                for kmode = 1:nmodes
                    Us{kmode} = ones(sizeX(kmode), lowerdims(kmode));
                end
        end
    end
else
    % initialisation as proposed in li14
    Us = cell(1, nmodes);
    for kmode = 1:nmodes
        Us{kmode} = ones(sizeX(kmode), lowerdims(kmode));
    end
end




%% run CMDA
errs = NaN(Tmax, 1);

if nargout>=4
    objfuncvals = NaN(nmodes*Tmax, 1);
end

if nargout>=5
    objfuncvals_traceratio = NaN(nmodes*Tmax, 1);
end

% calculate Xc - X for each class, where Xc is the class mean and X is the
% overall mean (stored in classmeandiffs) and Xcj - Xc where Xcj is the
% j'th observation from class c (stored in observationdiffs) and the number
% of observations from each class (stored in nis).
[classmeandiffs, observationdiffs, nis] = classbased_differences(Xs, classes);

[nclasses, ~, nmodes] = get_sizes(classmeandiffs, 1); % observations assumed to run along first mode

permute_vector = [2:(nmodes+1), 1]; % move observations to run along last mode
classmeandiffstensor = permute(classmeandiffs, permute_vector);
observationdiffstensor = permute(observationdiffs, permute_vector);

Rw = observationdiffstensor;
Rb = classwise_scalar_multiply(classmeandiffstensor, sqrt(nis));
% multiply all entries in classmeandiffstensor by the square root of the
% size of their class. When Rb is multiplied by its own transpose, the
% class sizes are automatically accounted for in the resulting sum.

stop = false;
iit = 0;
innerits = 0;
while ~stop && iit < Tmax
    iit = iit+1;
    oldUs = Us;
    for kmode = 1:nmodes
        innerits = innerits +1;
        
        permute_vector = 1:(nmodes+1);
        permute_vector(1) = kmode;
        permute_vector(kmode) = 1;
        
        othermodes = setdiff(1:nmodes, kmode);
        
        QtRb_mm = Rb;
        for othermode = othermodes
        QtRb_mm=tmult(QtRb_mm,Us{othermode}', othermode);
        end
        
        QtRb=reshape(permute(QtRb_mm, permute_vector),[sizeX(kmode),...
            prod(lowerdims(othermodes))*nclasses]);
        %QtRb=reshape(permute(QtRb_mm, [kmode, othermode, 3]),[sizeX(kmode),...
        %    lowerdims(othermode)*nclasses]);
        B = QtRb*QtRb';
        
        QtRw_mm = Rw;
        for othermode = othermodes
        QtRw_mm=tmult(QtRw_mm,Us{othermode}',othermode);
        end
        QtRw=reshape(permute(QtRw_mm, permute_vector),[sizeX(kmode),...
            prod(lowerdims(othermodes))*nobs]);
        W = QtRw*QtRw';
        
        
        [U, ~] = svd(W\B, 0);
        Us{kmode} = U(:, 1:lowerdims(kmode));
        
        if nargout >=4
            Btemp = Us{kmode}'*B*Us{kmode};
            Wtemp = Us{kmode}'*W*Us{kmode};
            objfuncvals(innerits) = -trace(Btemp)/trace(Wtemp);
        end
        if nargout >=5
            tempU.U1 = Us{1};
            tempU.U2 = Us{2};
            objfuncvals_traceratio(innerits) = tensorsldaobj_matrixdata(tempU,...
                classmeandiffs, observationdiffs, nis, lowerdims(1), lowerdims(2), ...
                Rw, Rb);
        end
    end
    
    % this is the stopping criterion proposed with CMDA in li14
    if usestoppingcrit && iit > 2
        errs(iit) = 0;
        for kmode = 1:nmodes
            errs(iit) = errs(iit) + norm(Us{kmode}*oldUs{kmode}' - eye(size(Us{kmode},1)), 'fro');
        end
        if errs(iit) <=tol
            stop = true;
        end
    end
    
end

if nargout >= 4
    objfuncvals = objfuncvals(1:(nmodes*iit));
end

if nargout >= 5
    objfuncvals_traceratio = objfuncvals_traceratio(1:(nmodes*iit));
end
if nargout >= 6
    Ys = tensor_projection(Xs, Us);
end

end