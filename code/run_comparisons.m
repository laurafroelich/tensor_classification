clear
addpath(genpath('code/'))
addpath('../matlab_additions/immoptibox/')

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
auc_lda = get_vectorised_lda_auc(x_train, x_test, y_train, y_test);
auc_dgtda = heuristic_project_predict(@DGTDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_dater = heuristic_project_predict(@DATER, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_datereig = heuristic_project_predict(@DATEReig, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_cmda = heuristic_project_predict(@CMDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_ManPDA = manifold_project_predict(@ManPDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_ManTDA = manifold_project_predict(@ManTDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_ManPDA_normsratio = manifold_project_predict(@ManPDA_normsratio, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_ManTDA_normsratio = manifold_project_predict(@ManTDA_normsratio, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_BDCA = direct_predict(@bilinear_logreg, x_train_cell, x_test_cell, y_train, y_test, ncomps);
auc_BDCA_tucker = direct_predict(@bilinear_logreg_tucker, x_train_cell, x_test_cell, y_train, y_test, ncomps);


