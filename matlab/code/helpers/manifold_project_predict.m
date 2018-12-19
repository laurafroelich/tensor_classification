function auc = manifold_project_predict(fun, x_train, x_test, ...
y_train, y_test, ncomps)
% auc = manifold_project_predict(fun, x_train, x_test, y_train, y_test, ncomps)
%
% Using one of the four projection methods optimised on the Stiefel
% manifold (@ManTDA_normsratio, @ManPDA_normsratio ', '@ManTDA, or
% @ManPDA), learn a projection, use this projection to project training and
% test data onto a lower dimensional representation. Then use these
% representations from the training data to train a multinomial regression
% classifier, and apply the classifier to the test data. Return the Area
% Under ROC Curve from the test data.
%
% Input:
%
% fun: manifold projection method (@ManTDA_normsratio, @ManPDA_normsratio ', 
%      '@ManTDA, or @ManPDA)
% x_train: cell array with one training observation in each cell.
% x_test: cell array with one test observation in each cell.
% y_train: vector with class labels. Only two classes currently supported,
%          these must be the labels 1 or 2.
% y_test: vector with class labels. Only two classes currently supported,
%          these must be the labels 1 or 2.
% ncomps: integer. Number of components in each mode.

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
predictions = project_and_predict(x_train, x_test, y_train, y_test, ncomps, Us);    
end

[~, ~, ~, auc] = perfcurve(y_test, predictions(:, 2), 2);



end