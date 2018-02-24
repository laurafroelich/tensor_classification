clear
addpath(genpath('code/'))

n_rows = 2;
n_cols = 3;
n_train_obs = 4;
n_test_obs = 10;

%x_train, y_train = 
[x_train, y_train] = simulate_data(n_train_obs, n_rows, n_cols);
[x_test, y_test] = simulate_data(n_test_obs, n_rows, n_cols);

x_train_mat = reshape(permute(x_train, [3, 1, 2]), [n_train_obs, n_rows*n_cols]);
x_test_mat = reshape(permute(x_test, [3, 1, 2]), [n_test_obs, n_rows*n_cols]);

lda_classifier = fitcdiscr(x_train_mat, y_train);
predictions = lda_classifier.predict(x_test_mat);


[~, ~, ~, auc] = perfcurve(y_test, predictions, 1);




