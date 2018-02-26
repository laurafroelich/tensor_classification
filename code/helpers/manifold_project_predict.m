function auc = manifold_project_predict(fun, x_train, x_test, ...
y_train, y_test, ncomps)
nrandinits = 5;
objfuncvals = NaN(1, nrandinits);
results = cell(1, nrandinits);

if isequal(fun, @ManTDA_normsratio) || isequal(fun, @ManPDA_normsratio) || ...
   isequal(fun, @ManTDA) || isequal(fun, @ManPDA)
else
        error(['manifold_project_predict.m: function handle given as ', ...
            'input must be one of @ManTDA_normsratio, @ManPDA_normsratio ', ...
            '@ManTDA, or @ManPDA.'])
end

varargs{1} = [ncomps, ncomps];

for iit=1:nrandinits
    [Us, outputs] = fun(x_train,...
        y_train, varargs{:});
    
    objfuncvals(iit) = outputs.fvals(end);
    
    results{iit}.Us = Us;
end

[~, best_solution] = min(objfuncvals);

Us = results{best_solution}.Us;

if isequal(fun, @ManPDA_normsratio) || isequal(fun, @ManPDA)
predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us, true);
else
predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Usxs);    
end

[~, ~, ~, auc] = perfcurve(y_test, predictions(:, 2), 2);



end