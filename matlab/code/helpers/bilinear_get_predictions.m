function predictions = bilinear_get_predictions(fun, x_train, x_test, ...
    y_train, ncomps, nrandinits)
% function predictions = bilinear_get_predictions(fun, x_train, x_test, ...
%    y_train, ncomps, nrandinits)
%
% Use a bilinear model to get predictions.
%
% This function uses ´x_train´ and ´y_train´ to optimise the parameters
% for a bilinear predictive model, which is then applied to ´x_test´ for 
% prediction.
%
% Input:
% 
% fun: function handle
%      The optimisation function to use, which determines whether a 
%      PARAFAC or a Tucker structure will be used in the model. Allowed
%      values are: 
%      - @bilinear_logreg, to use the PARAFAC structure
%      - @bilinear_logreg_tucker to use the Tucker structure
% x_train: 3D array or cell array containing equally sized matrices
%          If a cell array, each cell must contain a matrix corresponding
%          to an observation. The first dimension of the cell array must
%          be the number of observations, and the second must be 1.
%          If a 3D array, observations must run along the first mode, such
%          that the length of the first mode is the number of observations.
% x_test: cell array
%         Each cell must contain a matrix corresponding to an observation. 
%         The first dimension of the cell array must be the number of 
%         observations, and the second must be 1. If ´x_test´ is not a cell
%         array, ´x_train´ will be used in its place, and predictions will
%         be for the training data.
% y_train: vector
%          Class labels for each observation in the training data. The
%          class labels must be the integers 1 and 2.
% ncomps: integer
%         The number of components to use in the lower space. The same
%         number is used for each mode in the current implementation.
% nrandinits: integer
%             The number of random initialisations to start from. The
%             initialisation that resulted in the lowest objective function
%             value is used in the final result.
%
% See the following reference for details on the bilinear model with the
% PARAFAC structure:
% M. Dyrholm, C. Christoforou, and L. C. Parra, 
% ?Bilinear discriminant component analysis,? 
% The Journal of Machine Learning Research, vol. 8, pp. 1097?1111, 2007.
%
% See the following reference for a description of our modification of the
% bilinear model, with the Tucker structure:
% Frølich, L., Andersen, T. & Mørup, M. 
% Rigorous optimisation of multilinear discriminant analysis with Tucker 
% and PARAFAC structures. BMC Bioinformatics 19, 197 (2018) 
% doi:10.1186/s12859-018-2188-0 (https://rdcu.be/bWwIL)
    
if isequal(fun, @bilinear_logreg) || isequal(fun, @bilinear_logreg_tucker)
else
    error(['direct_predict.m: function handle given as ', ...
        'input must be one of @bilinear_logreg, or @bilinear_logreg_tucker.'])
end

nits=100;

if isa(x_train, 'cell')
    Xsmat_train = cell_array_to_nd_array(x_train);
else
    Xsmat_train = x_train;
end
permute_vector = [2, 3, 1];
Xsmat_train = permute(Xsmat_train, permute_vector);

if isa(x_test, 'cell')
    Xsmat_test = cell_array_to_nd_array(x_test);
else
    Xsmat_test = x_train;
end

Xsmat_test = permute(Xsmat_test, permute_vector);
results = cell(nrandinits, 1);
objfuncvals = NaN(1, nrandinits);
opts.maxeval = nits;
for iit = 1:nrandinits
    [I, J, ~] = size(Xsmat_train);
    if isequal(fun, @bilinear_logreg_tucker)
        nparams = 1+I*ncomps+J*ncomps+ncomps^2-ncomps;
    else
        nparams = 1+I*ncomps+J*ncomps;
    end
    x0 = randn(nparams, 1)*1e-5;
    
    %https://se.mathworks.com/matlabcentral/answers/149059-how-do-i-invoke-a-shadowed-core-matlab-function-not-built-in-from-an-overloaded-function-of-same
    this = ['linesearch', '.m']; % the name of function in MATLAB we are shadowing
    list = which(this, '-all'); % find all the functions which shadow it
    f = strfind(list, 'immoptibox'); % locate 1st in list under matlabroot
    correct_cell = find(~cellfun(@isempty,f));
    linesearch_path = list{correct_cell};
    [immoptibox_path, ~, ~] = fileparts(linesearch_path);
    
    here = cd(immoptibox_path); % temporarily switch to the containing folder
    tic
    [X, info, perf] = ucminf(fun, x0, opts, [], Xsmat_train, y_train-1, ncomps);
    t = toc;
    
    cd(here); % go back to where we came from
    
    finalxvector = X(:, end);
    
    if isequal(fun, @bilinear_logreg_tucker)
        w0 = finalxvector(1);
        Us{1} = reshape(finalxvector(2:(2+I*ncomps-1)),I,ncomps); % U in bilinear logreg code
        Us{2} = reshape(finalxvector((2+I*ncomps):(2+I*ncomps+J*ncomps-1)),...
            J,ncomps); % V in bilinear logreg code
        w = finalxvector(2+I*ncomps+J*ncomps:end);
        results{iit}.w = w;
    else
        w0 = finalxvector(1);
        Us{1} = reshape(finalxvector(2:(2+I*ncomps-1)),I,ncomps); % U in bilinear logreg code
        Us{2} = reshape(finalxvector((2+I*ncomps):end),J,ncomps); % V in bilinear logreg code
    end
    
    results{iit}.w0 = w0;
    results{iit}.Us = Us;
    results{iit}.X = X;
    results{iit}.info = info;
    results{iit}.perf = perf;
    results{iit}.t = t;
    if ~isempty(perf)
        results{iit}.objfuncvals = perf.f;
    end
    objfuncvals(iit) = info(1);
end

[~, bestit] = max(objfuncvals);



Us = results{bestit}.Us;
w0 = results{bestit}.w0;
XsU = tmult(Xsmat_test, Us{1}', 1);


if isequal(fun, @bilinear_logreg_tucker)
    W=nan(ncomps);
    ignore=logical(eye(ncomps));
    W(~ignore) = results{bestit}.w;
    W(ignore)=1;
    psiXs = squeeze(sum(sum(bsxfun(@times, tmult(XsU, Us{2}', 2), W),1),2));
    
else
    psiXs = squeeze(sum(sum(bsxfun(@times, XsU, Us{2}'),1),2)); % psiXs squeezed
end


linfunc = psiXs + w0; % vector of length ntrials
predprobsclasstwo = (1./(1+exp(-linfunc)));
predictions = [1-predprobsclasstwo predprobsclasstwo];
end