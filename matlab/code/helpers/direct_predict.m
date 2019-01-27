function auc = direct_predict(fun, x_train, x_test, ...
    y_train, y_test, ncomps)
if isequal(fun, @bilinear_logreg) || isequal(fun, @bilinear_logreg_tucker)
else
    error(['direct_predict.m: function handle given as ', ...
        'input must be one of @bilinear_logreg, or @bilinear_logreg_tucker.'])
end

nrandinits = 5;

predictions = bilinear_get_predictions(fun, x_train, x_test, ...
    y_train, ncomps, nrandinits);

[~, ~, ~, auc] = perfcurve(y_test, predictions(:, 2), 2);

%%%%%%%%









end