function [Us, objfuncval, Ys] = HODA(Xs, classes, varargin)
% [Us, objfuncval, Ys] = HODA(Xs, classes, varargin)
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
% page 15 in
% http://www.bsp.brain.riken.jp/publications/2010/IEICE_NOLTA_Phan-Cichocki-corr.pdf

if isa(Xs, 'cell')
    Xs = cell_array_to_nd_array(Xs);
end
    
sizeXs = size(Xs);
sizeX = sizeXs(2:end); % observations assumed to run along first mode
nmodes = length(sizeX);


if length(sizeX) > 2
    error(['HODA.m: Input data has more than two dimensions. '...
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

Us = cell(1, nmodes);
for kmode = 1:nmodes
    Us{kmode} = orth(randn(sizeX(kmode), lowerdims(kmode)));
end

% calculate Xc - X for each class, where Xc is the class mean and X is the
% overall mean (stored in classmeandiffs) and Xcj - Xc where Xcj is the
% j'th observation from class c (stored in observationdiffs) and the number
% of observations from each class (stored in nis).
[classmeandiffs, observationdiffs, nis] = classbased_differences(Xs, classes);


sizeobs = size(classmeandiffs);
I = sizeobs(2);
J = sizeobs(3);
nclasses = sizeobs(1);
nsamples = size(observationdiffs, 1);

permute_vector = [2:(length(sizeobs)), 1];
classmeandiffstensor = permute(classmeandiffs, permute_vector);
observationdiffstensor = permute(observationdiffs, permute_vector);
Xs = permute(Xs, permute_vector);

Rw = observationdiffstensor;
Rb = classmeandiffstensor.*permute(repmat(sqrt(nis), I,1,J), [1 3 2]);
% multiply all entries in classmeandiffstensor by the square root of the
% size of their class. When Rb is multiplied by its own transpose, the
% class sizes are automatically accounted for in the resulting sum.

maxits = 1000;


its = 0;
while true && its < maxits
    its = its + 1;
    oldUs = Us;
    difference = 0;
    for kmode = 1:nmodes
        othermode = setdiff(1:2, kmode);
        
        QtRb_mm=tmult(Rb,Us{othermode}', othermode);
        QtRb=reshape(permute(QtRb_mm, [kmode, othermode, 3]),[sizeX(kmode),...
            lowerdims(othermode)*nclasses]);
        B = QtRb*QtRb';
        
        QtRw_mm=tmult(Rw,Us{othermode}',othermode);
        QtRw=reshape(permute(QtRw_mm, [kmode, othermode, 3]),[sizeX(kmode),...
            lowerdims(othermode)*nsamples]);
        W = QtRw*QtRw';
        
        phi = trace(Us{kmode}' * B * Us{kmode})/trace(Us{kmode}' * W * Us{kmode});
        
        [U, ~] = eigs(B-phi*W, lowerdims(kmode));
        UUt = U*U';
        Xs_minus_n = matricizing(tmult(Xs, Us{othermode}', othermode), kmode);
        [Us{kmode}, ~] = eigs(UUt * (Xs_minus_n*Xs_minus_n') * UUt, lowerdims(kmode));
        
        difference = norm(Us{kmode}-oldUs{kmode}, 'fro')/(lowerdims(kmode)*sum(sizeX));
    end
    
    if difference<1e-6
        break
    end
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