function auc = get_vectorised_lda_auc(x_train, x_test, y_train, y_test)

x_size = size(x_train);
n_rows = x_size(1);
n_cols = x_size(2);
n_train_obs = x_size(end);
y_test_size = size(y_test);
n_test_obs = y_test_size(1);

x_train_mat = reshape(permute(x_train, [3, 1, 2]), [n_train_obs, n_rows*n_cols]);
x_test_mat = reshape(permute(x_test, [3, 1, 2]), [n_test_obs, n_rows*n_cols]);

lda_classifier = fitcdiscr(x_train_mat, y_train);
predictions = lda_classifier.predict(x_test_mat);


[~, ~, ~, auc] = perfcurve(y_test, predictions, 2);
end