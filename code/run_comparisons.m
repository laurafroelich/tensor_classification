clear
addpath(genpath('code/'))

n_rows = 2;
n_cols = 3;
n_train_obs = 5;
n_test_obs = 10;

%x_train, y_train = 
[x_train, y_train] = simulate_data(n_train_obs, n_rows, n_cols);
x_train_cell = mat_to_cell(x_train);

[x_test, y_test] = simulate_data(n_test_obs, n_rows, n_cols);
x_test_cell = mat_to_cell(x_test);

ncomps = 2;
auc_lda = get_vectorised_lda_auc(x_train, x_test, y_train, y_test)
auc_dgtda = get_dgtda_auc(x_train_cell, x_test_cell, y_train, y_test, ncomps)

