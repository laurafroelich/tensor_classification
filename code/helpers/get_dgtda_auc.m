function auc = get_dgtda_auc(x_train, x_test, y_train, y_test, ncomps)
nrandinits = 5;
objfuncvals = NaN(1, nrandinits);
results = cell(1, nrandinits);
for iit=1:nrandinits
    [Us, objfuncval] = DGTDA(x_train,...
        y_train, [ncomps ncomps ncomps]);
    
    objfuncvals(iit) = objfuncval;
    
    results{iit}.Us = Us;
end

[~, best_solution] = min(objfuncvals);

Us = results{best_solution}.Us;

predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us);

[~, ~, ~, auc] = perfcurve(y_test, predictions(:, 2), 2);



end