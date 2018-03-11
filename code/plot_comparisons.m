addpath(genpath('code/'))
addpath('../matlab_additions/immoptibox/')
run('colsandlinestyles.m')
savefig = false;

load('results/aucs.mat')
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
    'Location', 'SouthWest')
title('AUC vs number of training samples', 'FontSize', titlefontsize)
set(gca, 'FontSize', gcafontsize)
xlabel('Number of training samples')
ylabel('Area Under ROC Curve')
set(gcf,'color','w');

if savefig
    export_fig(figh, 'figures/ManTDA_sr_vs_nonMan_structured_noise_high_aucs_little_training.pdf', '-pdf')
end
