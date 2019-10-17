function [Us, iit, errs, objfuncvals, objfuncvals_traceratio, Ys] = DATEReig(Xs, classes, varargin)
% [Us, iit, errs, objfuncvals, Ys] = DATER(Xs, classes, varargin)
% Mandatory input:
% Xs:           Cell array containing the observed tensors.
% classes:      Vector containing class labels. Classes must be sequential
%               numbers starting from one.
%
% Optional input:
% varargin{1}:  Scalar giving the maximal number of iterations.
%               Default: 100.
% varargin{2}:  Vector containing the number of components to estimate for
%               each mode. Default: size of observations.
% varargin{3}:  Boolean indicating whether or not convergence criteria
%               should be allowed to terminate iterations.
%               Default: true.
% varargin{4}:  Cell array of initial projection matrices or the string
%               'randinit' indicating initialisation with random
%               orthogonal matrices. Default: identity matrices as
%               proposed in yan05.
%
% Output:
% Us:           Cell array containing the projection matrices found.
% iit:          The number of outer iterations performed (an outer 
%               iteration consists of one optimisation of each mode).
% errs:         The values from the stopping criterion proposed in yan05, 
%               one for each outer iteration.
% objfuncvals:  Values of the objective function that DATER tries to
%               optimise. Values are given for all inner iterations. An
%               inner iteration consists of the optimisation of one mode.
% objfuncvals_traceratio:   Values of the objective function
%                           tr((U'WU)^(-1)U'BU). W is
%                           the within-class scatter matrix while B is the
%                           between-class scatter matrix.
% Ys:           Projections of the original data in Xs projected onto the
%               final projection matrices in Us.
%
% vis05:
%   Visani, Muriel and Garcia, Christophe and Jolion, J-M
%   'Normalized radial basis function networks and bilinear discriminant analysis for face recognition'
%   Advanced Video and Signal Based Surveillance, 2005. AVSS 2005. IEEE Conference on

% yan05:
%   S. Yan, D. Xu, Q. Yang, L. Zhang, X. Tang, and H. Zhang
%   'Discriminant Analysis with Tensor Representation'
%   IEEE Int. Conf. on Computer Vision and Pattern Recognition, 2005

%% read input and set parameters
if isa(Xs, 'cell')
    Xs = cell_array_to_nd_array(Xs);
end
    
[nobs, sizeX, nmodes] = get_sizes(Xs, 1); % observations assumed to run along first mode

tol=1e-6;

if length(varargin)>=4 && ~isempty(varargin{4})
    Us = varargin{4};
    if ischar(Us)
        if ~(strcmp(Us, 'randinit') || strcmp(Us, 'ones') || ...
                strcmp(Us, 'identity'))
                warning(['DATEReig.m: initialisation method not recognised, '...
                    'initialising with identity matrices as proposed in yan05 (see help for citation)'])
                % initialisation as proposed in yan05
                varargin{4} = 'identity';
        end
    end
else
    varargin{4} = 'identity';
end

[Tmax, lowerdims, usestoppingcrit, Us] = parse_varargin(sizeX, nmodes, varargin{:});


%% run DATEReig vis05
errs = NaN(Tmax, 1);
if nargout >=4
    objfuncvals = NaN(nmodes*Tmax, 1); % 
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

innerits = 0;
for iit = 1:Tmax
    
    oldUs = Us;
    for kmode = 1:nmodes
        innerits = innerits +1;
        between_class_scatter = get_mode_specific_scatter_matrix(Rb, kmode, lowerdims, Us, ...
            nclasses, nmodes, sizeX);
        
        within_class_scatter = get_mode_specific_scatter_matrix(Rw, kmode, lowerdims, Us, ...
            nobs, nmodes, sizeX);
         
        [U, eigvals] = eig(within_class_scatter\between_class_scatter);
        
        eigvals = diag(eigvals);
        [~, sortedinds] = sort(eigvals, 'descend');
        
        Us{kmode} = U(:, sortedinds(1:lowerdims(kmode)));
        
        if nargout >=4
            Btemp = Us{kmode}'*B*Us{kmode};
            Wtemp = Us{kmode}'*W*Us{kmode};
            objfuncvals(innerits) = -trace(Btemp)/trace(Wtemp);
        end
        if nargout >=5
            tempU.U1 = Us{1};
            tempU.U2 = Us{2};
            objfuncvals_traceratio(innerits) = tensorsldaobj_matrixdata(tempU,...
                classmeandiffs, observationdiffs, nis, lowerdims, ...
                Rw, Rb);
        end
    end
    
    % stopping criterion proposed with DATER in yan05
    if usestoppingcrit && (iit > 2)
        stopnow = true;
        for kmode = 1:nmodes
            errs(iit) = norm(Us{kmode} - oldUs{kmode});
            stopnow = stopnow && (norm(Us{kmode} -...
                oldUs{kmode})<sizeX(kmode)*lowerdims(kmode)*tol);
        end
        if stopnow
            break
        end
    end
end

if nargout >=4
    objfuncvals = objfuncvals(1:(nmodes*iit));
end

if nargout >= 5
    objfuncvals_traceratio = objfuncvals_traceratio(1:(nmodes*iit));
end
if nargout >= 6
    Ys = tensor_projection(Xs, Us);
end

end