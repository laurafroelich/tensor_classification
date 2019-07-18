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
    
[nobs, sizeX, nmodes] = get_sizes(Xs, 1); % observations assumed to run along first mode

[~, lowerdims, ~, Us] = parse_varargin(sizeX, nmodes, [], varargin{1}, ...
    [], 'orth');

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

maxits = 1000;

its = 0;
while true && its < maxits
    its = its + 1;
    oldUs = Us;
    difference = 0;
    for kmode = 1:nmodes
        
        between_class_scatter = get_mode_specific_scatter_matrix(Rb, kmode, lowerdims, Us, ...
            nclasses, nmodes, sizeX);
        
        within_class_scatter = get_mode_specific_scatter_matrix(Rw, kmode, lowerdims, Us, ...
            nobs, nmodes, sizeX);
        
        
        phi = trace(Us{kmode}' * between_class_scatter * Us{kmode})/...
            trace(Us{kmode}' * within_class_scatter * Us{kmode});
        
        [U, ~] = eigs(between_class_scatter-phi*within_class_scatter, lowerdims(kmode));
        UUt = U*U';
        Xs_minus_n = Xs;
        for mode = setdiff(1:nmodes, kmode)
            Xs_minus_n = tmult(Xs_minus_n, Us{mode}', mode+1);
        end
        Xs_minus_n = matricizing(Xs_minus_n, kmode+1);
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