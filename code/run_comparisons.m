addpath(genpath('code/'))
addpath('../matlab_additions/immoptibox/')
run('colsandlinestyles.m')
savefig = false;

n_rows = 10;
n_cols = 80;
trainobs = 10*2.^(0:8); 
configs = {'high', 'medium', 'low'};
n_test_obs = 500;
n_true_comps = 3;
ncomps = 3;

for iconfig = 1:1
    [sigma, corestd, coreNoiseStd, UNoiseStd] = get_config(configs{iconfig});
    
    for i = 8:length(trainobs)
        for iit = 1:10
            n_train_obs = trainobs(i);
            
            [x, y] = simulate_data(n_train_obs+n_test_obs, n_rows, n_cols, ...
                n_true_comps, sigma, corestd, coreNoiseStd, UNoiseStd);
            
            c = cvpartition(y, 'HoldOut', n_test_obs/(n_train_obs+n_test_obs));
            x_train = x(:,:,c.training);
            y_train = y(c.training);
            x_test = x(:,:,c.test);
            y_test= y(c.test);
            
            x_train_cell = mat_to_cell(x_train);
            x_test_cell = mat_to_cell(x_test);
            
            
            auc_lda(iit, i, iconfig) = get_vectorised_lda_auc(x_train, x_test, y_train, y_test);
            auc_dgtda(iit, i, iconfig) = heuristic_project_predict(@DGTDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_dater(iit, i, iconfig) = heuristic_project_predict(@DATER, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_datereig(iit, i, iconfig) = heuristic_project_predict(@DATEReig, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_cmda(iit, i, iconfig) = heuristic_project_predict(@CMDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_ManPDA(iit, i, iconfig) = manifold_project_predict(@ManPDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_ManTDA(iit, i, iconfig) = manifold_project_predict(@ManTDA, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_ManPDA_normsratio(iit, i, iconfig) = manifold_project_predict(@ManPDA_normsratio, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_ManTDA_normsratio(iit, i, iconfig) = manifold_project_predict(@ManTDA_normsratio, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_BDCA(iit, i, iconfig) = direct_predict(@bilinear_logreg, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_BDCA_tucker(iit, i, iconfig) = direct_predict(@bilinear_logreg_tucker, x_train_cell, x_test_cell, y_train, y_test, ncomps);
            auc_tucker(iit, i, iconfig) = decompose_predict(@tucker_decompose, x_train, x_test, y_train, y_test, ncomps);
            auc_parafac(iit, i, iconfig) = decompose_predict(@parafac_decompose, x_train, x_test, y_train, y_test, ncomps);
            auc_tucker2(iit, i, iconfig) = decompose_predict(@tucker2_decompose, x_train, x_test, y_train, y_test, ncomps);
            %auc_parafac2(i) = decompose_predict(@parafac2_decompose, x_train, x_test, y_train, y_test, ncomps);
            save(['results/aucs_iit_', num2str(iit), '.mat'])
        end
    end
end
save('results/aucs.mat')
