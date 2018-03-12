function auc = heuristic_project_predict(fun, x_train, x_test, y_train, y_test, ncomps)
nrandinits = 5;
objfuncvals = NaN(1, nrandinits);
results = cell(1, nrandinits);

if isequal(fun, @DGTDA) || isequal(fun, @HODA)
    varargs{1} = [ncomps, ncomps];
else
    if isequal(fun, @CMDA) || isequal(fun, @DATER) || ...
        isequal(fun, @DATEReig) || isequal(fun, @HODA)
        varargs{2} = [ncomps, ncomps];
        varargs{4} = 'randinit';
    else
        error(['heuristic_project_predict.m: function handle given as ', ...
            'input must be one of @DGTDA, @DATER, @DATEReig, @HODA, or @CMDA.'])
    end
end

for iit=1:nrandinits
    [Us, objfuncval] = fun(x_train,...
        y_train, varargs{:});
    
    objfuncvals(iit) = objfuncval;
    
    results{iit}.Us = Us;
end

[~, best_solution] = min(objfuncvals);

Us = results{best_solution}.Us;

predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us);

[~, ~, ~, auc] = perfcurve(y_test, predictions(:, 2), 2);



end