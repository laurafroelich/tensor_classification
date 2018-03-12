addpath(genpath('code/'))
addpath('../matlab_additions/immoptibox/')
run('colsandlinestyles.m')
savefig = false;

load('results/aucs.mat')
%% log x axis plots of mean
for iconfig = 1:3

figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
semilogx(trainobs, mean(auc_lda(:,:,iconfig), 1))
hold on;
semilogx(trainobs, mean(auc_dgtda(:,:,iconfig), 1), 'o-')
semilogx(trainobs, mean(auc_dater(:,:,iconfig), 1))
semilogx(trainobs, mean(auc_datereig(:,:,iconfig), 1))
semilogx(trainobs, mean(auc_cmda(:,:,iconfig), 1))
semilogx(trainobs, mean(auc_hoda(:,:,iconfig), 1), 'o-')


semilogx(trainobs, mean(auc_ManPDA(:,:,iconfig), 1), 'r-', 'LineWidth', 2)
semilogx(trainobs, mean(auc_ManTDA(:,:,iconfig), 1), 'g-', 'LineWidth', 2)
semilogx(trainobs, mean(auc_ManPDA_normsratio(:,:,iconfig), 1), 'b-', 'LineWidth', 2)
semilogx(trainobs, mean(auc_ManTDA_normsratio(:,:,iconfig), 1), 'k-', 'LineWidth', 2)
semilogx(trainobs, mean(auc_tucker(:,:,iconfig), 1), '-x')
semilogx(trainobs, mean(auc_parafac(:,:,iconfig), 1), '-x')
semilogx(trainobs, mean(auc_tucker2(:,:,iconfig), 1), '-x')
%plot(trainobs, auc_parafac2, '-x')
legend('LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA',...
    'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
    'Tucker', 'PARAFAC', 'Tucker2',... % 'parafac2', ...
    'Location', 'NorthEast')
title('Mean AUC vs number of training samples', 'FontSize', titlefontsize)
set(gca, 'FontSize', gcafontsize)
xlabel('Number of training samples')
ylabel('Mean Area Under ROC Curve')
set(gcf,'color','w');
ylim([0.4, 1])
if savefig
    export_fig(figh, ['figures/auc_avgss_signal_strengh_', configs{iconfig}, '.pdf'], '-pdf')
end
end


%% log x axis plots of std

for iconfig = 1:3

figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
semilogx(trainobs, std(auc_lda(:,:,iconfig), [], 1))
hold on;
semilogx(trainobs, std(auc_dgtda(:,:,iconfig), [], 1))
semilogx(trainobs, std(auc_dater(:,:,iconfig), [], 1))
semilogx(trainobs, std(auc_datereig(:,:,iconfig), [], 1))
semilogx(trainobs, std(auc_cmda(:,:,iconfig), [], 1))


semilogx(trainobs, std(auc_ManPDA(:,:,iconfig), [], 1), 'r-', 'LineWidth', 2)
semilogx(trainobs, std(auc_ManTDA(:,:,iconfig), [], 1), 'g-', 'LineWidth', 2)
semilogx(trainobs, std(auc_ManPDA_normsratio(:,:,iconfig), [], 1), 'b-', 'LineWidth', 2)
semilogx(trainobs, std(auc_ManTDA_normsratio(:,:,iconfig), [], 1), 'k-', 'LineWidth', 2)
semilogx(trainobs, std(auc_tucker(:,:,iconfig), [], 1), '-x')
semilogx(trainobs, std(auc_parafac(:,:,iconfig), [], 1), '-x')
semilogx(trainobs, std(auc_tucker2(:,:,iconfig), [], 1), '-x')
%plot(trainobs, auc_parafac2, '-x')
legend('LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA',...
    'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
    'Tucker', 'PARAFAC', 'Tucker2',... % 'parafac2', ...
    'Location', 'NorthEast')
title('Standard deviation of AUC vs number of training samples', 'FontSize', titlefontsize)
set(gca, 'FontSize', gcafontsize)
xlabel('Number of training samples')
ylabel('Standard deviation of Area Under ROC Curve')
set(gcf,'color','w');
ylim([0, 0.3])

if savefig
    export_fig(figh, ['figures/auc_stds_signal_strengh_', configs{iconfig}, '.pdf'], '-pdf')
end
end