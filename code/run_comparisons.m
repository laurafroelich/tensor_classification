clear
addpath(genpath('code/'))
addpath('../matlab_additions/immoptibox/')
run('colsandlinestyles.m')

n_rows = 10;
n_cols = 80;

savefig = false;

trainobs = 100:200:2000; %[50:50:250]; %%[50, 100, 200, 400, 800, 1600];% [200, 250, 300];%, 100, 150, 200, 250, 300, 400, 500];%, 75, 100];%, 150, 200, 250, 300, 350, 400, 500, 600];

for i = 1:length(trainobs)
    ncomps = 3;
    
    n_train_obs = trainobs(i);
    n_test_obs = 500;
    n_true_comps = 3;
    
    rng(0)
    [x, y] = simulate_data(n_train_obs+n_test_obs, n_rows, n_cols, n_true_comps);
    
    c = cvpartition(y, 'HoldOut', n_test_obs/(n_train_obs+n_test_obs));
    x_train = x(:,:,c.training);
    y_train = y(c.training);
    x_test = x(:,:,c.test);
    y_test= y(c.test);
    
    x_train_cell = mat_to_cell(x_train);
    x_test_cell = mat_to_cell(x_test);
    
    
    auc_lda(i) = get_vectorised_lda_auc(x_train, x_test, y_train, y_test);
    auc_dgtda(i) = heuristic_project_predict(@DGTDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_dater(i) = heuristic_project_predict(@DATER, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_datereig(i) = heuristic_project_predict(@DATEReig, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_cmda(i) = heuristic_project_predict(@CMDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_ManPDA(i) = manifold_project_predict(@ManPDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_ManTDA(i) = manifold_project_predict(@ManTDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_ManPDA_normsratio(i) = manifold_project_predict(@ManPDA_normsratio, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_ManTDA_normsratio(i) = manifold_project_predict(@ManTDA_normsratio, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_BDCA(i) = direct_predict(@bilinear_logreg, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_BDCA_tucker(i) = direct_predict(@bilinear_logreg_tucker, x_train_cell, x_test_cell, y_train, y_test, ncomps);
    auc_tucker(i) = decompose_predict(@tucker_decompose, x_train, x_test, y_train, y_test, ncomps);
    auc_parafac(i) = decompose_predict(@parafac_decompose, x_train, x_test, y_train, y_test, ncomps);
    auc_tucker2(i) = decompose_predict(@tucker2_decompose, x_train, x_test, y_train, y_test, ncomps);
    %auc_parafac2(i) = decompose_predict(@parafac2_decompose, x_train, x_test, y_train, y_test, ncomps);
end
%% linear x axis plots
close all

plot(trainobs, auc_lda)
hold on;
plot(trainobs, auc_dgtda)
plot(trainobs, auc_dater)
plot(trainobs, auc_datereig)
plot(trainobs, auc_cmda)
legend('lda', 'dgtda', 'dater', 'datereig', 'cmda')
ylim([0, 1])

figure()
plot(trainobs, auc_ManPDA, '-x')
hold on;
plot(trainobs, auc_ManTDA, '-x')
plot(trainobs, auc_ManPDA_normsratio, '-x')
plot(trainobs, auc_ManTDA_normsratio, '-x')
plot(trainobs, auc_BDCA, '-x')
plot(trainobs, auc_BDCA_tucker, '-x')
legend('ManPDA', 'ManTDA', 'ManPDA\_normsratio', ...
    'ManTDA\_normsratio', 'BDCA', 'BDCA\_tucker')
ylim([0, 1])


figure()
plot(trainobs, auc_tucker, '-x')
hold on;
plot(trainobs, auc_parafac, '-x')
plot(trainobs, auc_tucker2, '-x')
%plot(trainobs, auc_parafac2, '-x')
legend('tucker', 'parafac', 'tucker2')%, 'parafac2')
ylim([0, 1])
%legend('lda', 'dgtda', 'dater', 'datereig', 'cmda', 'ManPDA', 'ManTDA', 'ManPDA_normsratio', ...
%    'ManTDA_normsratio', 'BDCA', 'BDCA_tucker', 'tucker', 'parafac', 'tucker2', 'parafac2')

%% linear x axis plots 


figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(trainobs, auc_lda)
hold on;
plot(trainobs, auc_dgtda)
plot(trainobs, auc_dater)
plot(trainobs, auc_datereig)
plot(trainobs, auc_cmda)

plot(trainobs, auc_ManTDA_normsratio, 'k-', 'LineWidth', 2)
plot(trainobs, auc_tucker, '-x')
plot(trainobs, auc_parafac, '-x')
plot(trainobs, auc_tucker2, '-x')
%plot(trainobs, auc_parafac2, '-x')
legend('LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA',...
    'ManTDA\_sr', 'Tucker', 'PARAFAC', 'Tucker2',... % 'parafac2', ...
    'Location', 'NorthWest')
title('AUC vs number of training samples', 'FontSize', titlefontsize)
set(gca, 'FontSize', gcafontsize)
xlabel('Number of training samples')
ylabel('Area Under ROC Curve')
set(gcf,'color','w');
    
if savefig
    export_fig(figh, 'figures/ManTDA_sr_vs_nonMan_structured_noise_hihg_aucs.pdf', '-pdf')
end