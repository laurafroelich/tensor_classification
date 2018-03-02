function auc = decompose_predict(fun, x_train, x_test, y_train, y_test, ncomps)
nrandinits = 5;
objfuncvals = NaN(1, nrandinits);
results = cell(1, nrandinits);

if isequal(fun, @tucker_decompose) || isequal(fun, @tucker2_decompose) || ...
        isequal(fun, @parafac_decompose) || isequal(fun, @parafac2_decompose)
else
    error(['decompose_predict.m: function handle given as ', ...
        'input must be one of @tucker_decompose, @tucker2_decompose ', ...
        '@parafac_decompose, or @parafac2_decompose.'])
end

X = cat(3, x_train, x_test);

for iit=1:nrandinits
    [As, objfuncvals(iit)] = fun(X, ncomps);
    
    results{iit}.As = As;
end

[~, best_solution] = max(objfuncvals);

if isequal(fun, @tucker_decompose) || isequal(fun, @parafac_decompose) || ...
    isequal(fun, @parafac2_decompose)
    trial_strengths = results{best_solution}.As{3};
    n_train = length(y_train);
    training_data = trial_strengths(1:n_train, :);
    test_data = trial_strengths((n_train+1):end, :);
    
    mnr_mod = mnrfit(training_data, y_train);
    predictions = mnrval(mnr_mod, test_data);
end

if isequal(fun, @tucker2_decompose)
    Us = {results{best_solution}.As{1:2}};
    predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us);
    
end

[~, ~, ~, auc] = perfcurve(y_test, predictions(:, 2), 2);



end