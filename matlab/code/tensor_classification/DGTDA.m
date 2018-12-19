function [Us, objfuncval, Ys] = DGTDA(Xs, classes, varargin)
% [Us, objfuncval, Ys] = DGTDA(Xs, classes, varargin)
% Mandatory input:
% Xs:           Cell array containing the observed tensors.
% classes:      Vector containing class labels. Classes must be sequential
%               numbers starting from one.
%
% Optional input:
% varargin{1}:  Vector containing the number of components to estimate for
%               each mode. Default: size of observations.
%
% Output:
% Us:           Cell array containing the projection matrices found.
% objfuncvals:  Value of the objective function that DGTDA tries to
%               optimise, with the final projection matrices in Us.
%
% li14:
%   Q. Li and D. Schonfeld,
%   'Multilinear Discriminant Analysis for Higher-Order Tensor
%   Data Classification',
%   IEEE Transactions on Pattern Analysis and Machine Intelligence


Xsample1 = Xs{1};
sizeX = size(Xsample1);
nmodes = length(sizeX);
nsamples = length(Xs);
nclasses = length(unique(classes));

if length(sizeX) > 2
    error(['DGTDA.m: Input data has more than two dimensions. '...
        'This function is only customised for two-dimensional (i.e. matrix) data.'])
end

if isempty(varargin) || isempty(varargin{1})
    lowerdims = sizeX;
else
    lowerdims = varargin{1};
end

if isempty(lowerdims)
    lowerdims = sizeX;
end

% calculate Xc - X for each class, where Xc is the class mean and X is the
% overall mean (stored in classmeandiffs) and Xcj - Xc where Xcj is the
% j'th observation from class c (stored in observationdiffs) and the number
% of observations from each class (stored in nis).
[classmeandiffs, observationdiffs, nis] = classbased_differences(Xs, classes);

Bs = cell(nmodes, 1);
for nmode = 1:nmodes
    diffmatricised = matricizing(classmeandiffs{1}, nmode);
    Bs{nmode} = nis(1)*(diffmatricised*diffmatricised');
    for iclass=1:nclasses
        diffmatricised = matricizing(classmeandiffs{iclass}, nmode);
        Bs{nmode} = Bs{nmode} + nis(iclass)*(diffmatricised*diffmatricised');
    end
end


Ws = cell(nmodes, 1);
for nmode = 1:nmodes
    Ws{nmode} = zeros(size(Bs{nmode}));
    for isample=1:nsamples
        diffmatricised = matricizing(observationdiffs{isample}, nmode);
        Ws{nmode} = Ws{nmode} + (diffmatricised*diffmatricised');
    end
end

Us = cell(nmodes, 1);
for nmode = 1:nmodes
    %eigvals = eig(Ws{nmode}\Bs{nmode});
    singvals = svd(Ws{nmode}\Bs{nmode});
    eta = max(singvals);
    [U, ~, ~] = svd(Bs{nmode} - eta*Ws{nmode});
    Us{nmode} = U(:, 1:lowerdims(nmode));
end

if nargout >= 2
    
    Ustruct.U1 = Us{1};
    Ustruct.U2 = Us{2};
    objfuncval = tensorsldaobj_matrixdata_normsratio(Ustruct, classmeandiffs, observationdiffs, nis, lowerdims(1), lowerdims(2));
    if nargout >= 3
        Ys = tensor_projection(Xs, Us);
    end
end

end